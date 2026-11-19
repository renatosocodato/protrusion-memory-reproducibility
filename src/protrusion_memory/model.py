"""Exact linear paired model and its calibrated reporter comparison."""

from dataclasses import dataclass
import math
from numbers import Real

import numpy as np
from scipy.linalg import expm
from scipy.special import ndtr

CRITICAL = 1.959963984540054


def finite_number(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")


def integer(value, name, least):
    if type(value) is not int or value < least:
        raise ValueError(f"{name} must be an integer >= {least}")


@dataclass(frozen=True)
class Pulse:
    start_s: int
    stop_s: int
    amplitude: float

    def __post_init__(self):
        integer(self.start_s, "pulse start", 0)
        integer(self.stop_s, "pulse stop", self.start_s + 1)
        finite_number(self.amplitude, "pulse amplitude")


@dataclass(frozen=True)
class Model:
    fast_s: float = 8
    slow_s: float = 25
    deformation_s: float = 20
    duration_s: int = 300
    gate: Pulse = Pulse(20, 80, 0.3)
    actin: Pulse = Pulse(140, 200, 0.2)

    def __post_init__(self):
        for name in ("fast_s", "slow_s", "deformation_s"):
            finite_number(getattr(self, name), name)
        if not 0 < self.fast_s < self.slow_s or self.deformation_s <= 0:
            raise ValueError("require 0 < fast_s < slow_s and deformation_s > 0")
        integer(self.duration_s, "duration_s", 1)
        if self.duration_s < max(self.gate.stop_s, self.actin.stop_s):
            raise ValueError("pulse stop times must be within duration_s")

    @property
    def rho(self):
        return self.fast_s / self.slow_s

    @property
    def touching_bound(self):
        return (self.slow_s - self.fast_s) / self.slow_s


@dataclass(frozen=True)
class Reporter:
    response_s: float = 10
    cadence_s: int = 15
    noise_sd: float = 0.03
    mixing_bound: float = 0.20
    both_orders: bool = True

    def __post_init__(self):
        for name in ("response_s", "noise_sd", "mixing_bound"):
            finite_number(getattr(self, name), name)
        if self.response_s < 0 or self.noise_sd <= 0 or not 0 <= self.mixing_bound <= 1:
            raise ValueError("require response_s >= 0, noise_sd > 0, and 0 <= mixing_bound <= 1")
        integer(self.cadence_s, "cadence_s", 1)
        if type(self.both_orders) is not bool:
            raise ValueError("both_orders must be true or false")


def _inputs(model, order):
    if order not in ("uv", "vu"):
        raise ValueError("order must be uv (configured windows) or vu (swapped windows)")
    gate, actin = model.gate, model.actin
    gate_window, actin_window = (gate, actin) if order == "uv" else (actin, gate)
    times = np.arange(model.duration_s)
    return np.column_stack((
        gate.amplitude * ((gate_window.start_s <= times) & (times < gate_window.stop_s)),
        actin.amplitude * ((actin_window.start_s <= times) & (times < actin_window.stop_s)),
    ))


def _propagate(generator, inputs):
    # Augment the dynamic state with the constant input on each one-second interval.
    transition = expm(generator)
    size = generator.shape[0] - 2
    states = np.zeros((len(inputs) + 1, size))
    for index, pulse in enumerate(inputs):
        states[index + 1] = (transition @ np.r_[states[index], pulse])[:size]
    if not np.isfinite(states).all():
        raise ValueError("states are outside the numerical range")
    return states


def _balances(model, family):
    fast = model.fast_s
    slow = model.slow_s
    rho = model.rho
    generator = np.zeros((5, 5))
    if family == "A":
        generator[0, [0, 3]] = [-1 / fast, 1 / fast]
        generator[1, [0, 1, 4]] = [1 / slow, -1 / slow, 1 / slow]
    else:
        generator[0, [0, 3, 4]] = [-1 / slow, 1 / slow, model.touching_bound / slow]
        generator[1, [0, 1, 4]] = [1 / fast, -1 / fast, rho / fast]
    return generator


def simulate_pair(model=Model(), order="uv"):
    """Return A and B states (gate, actin, deformation) and the optical precursor in A."""
    inputs = _inputs(model, order)
    states = []
    for family in ("A", "B"):
        generator = _balances(model, family)
        generator[2, 1:3] = [1 / model.deformation_s, -1 / model.deformation_s]
        states.append(_propagate(generator, inputs))
    a, b = states
    optical_a = model.rho * a[:, 0] + model.touching_bound * a[:, 1]
    return a, b, optical_a


def filtered_contrast(model, response_s, order="uv"):
    """Propagate the common reporter template for the actin-subtracted signal."""
    finite_number(response_s, "response_s")
    if response_s < 0:
        raise ValueError("response_s must be >= 0")
    generator = _balances(model, "A")
    if response_s:
        generator[2, :3] = np.array([1, -1, -1]) / response_s
    values = _propagate(generator, _inputs(model, order))
    return values[:, 2] if response_s else values[:, 0] - values[:, 1]


def probabilities(gap_sigma):
    """Closest-mean Gaussian probabilities on a defined one-dimensional template."""
    threshold = max(gap_sigma / 2, CRITICAL)
    correct = float(ndtr(gap_sigma - threshold))
    wrong = float(ndtr(-threshold))
    inconclusive = 0.0 if threshold == gap_sigma / 2 else max(0.0, 1.0 - correct - wrong)
    return correct, wrong, inconclusive


def compare_reporter(model=Model(), reporter=Reporter()):
    """Profile worst-case separation over integer sampling phases and requested orders."""
    if reporter.cadence_s > model.duration_s + 1:
        raise ValueError("cadence_s must leave at least one observation at every phase")
    orders = ("uv", "vu") if reporter.both_orders else ("uv",)
    norms = []
    for order in orders:
        contrast = filtered_contrast(model, reporter.response_s, order)
        for phase in range(reporter.cadence_s):
            template = contrast[phase::reporter.cadence_s]
            norm = math.hypot(*template)
            if not math.isfinite(norm):
                raise ValueError("sampled template norm is outside the numerical range")
            norms.append((norm, order, phase))
    minimum, order, phase = min(norms)
    bound = reporter.mixing_bound
    if minimum == 0:
        status = "zero_template"
    elif abs(bound - model.touching_bound) <= 2 * math.ulp(model.touching_bound):
        status = "touching_means"
    elif bound > model.touching_bound:
        status = "overlapping_means"
    else:
        status = "separated_means"
    gap = max(0.0, (model.touching_bound - bound) * (minimum / reporter.noise_sd)) if (
        status == "separated_means"
    ) else 0.0
    if not math.isfinite(gap):
        raise ValueError("standardized separation is outside the numerical range")
    conditional = None
    null_declaration = None
    if status == "separated_means":
        correct, wrong, unresolved = probabilities(gap)
        conditional = {"correct": correct, "wrong": wrong, "inconclusive": unresolved}
    elif status == "touching_means":
        null_declaration = 2 * probabilities(0.0)[1]
    return {
        "status": status,
        "rho": model.rho,
        "touching_bound": model.touching_bound,
        "coefficient_intervals": {"A": [1 - bound, 1], "B": [model.rho * (1 - bound), model.rho]},
        "gap_sigma": gap,
        "least_informative": {"order": order, "phase_s": phase, "template_norm": minimum},
        "maximum_template_norm": max(norms)[0],
        "closest_mean_probabilities": conditional,
        "either_declaration_at_touching_mean": null_declaration,
    }
