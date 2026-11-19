"""Declarative deterministic ODEs, parameter interventions and observation maps."""

import ast
from bisect import bisect_right
from dataclasses import dataclass
from decimal import Decimal
import keyword
import math
import operator
import re

import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import quad, solve_ivp

SCHEMA = "ode-v1"
ACQUISITION_SCHEMA = "ode-v2"
SCHEMAS = (SCHEMA, ACQUISITION_SCHEMA)
RTOL, ATOL = 1e-9, 1e-11
BINARY = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
          ast.Div: operator.truediv, ast.Pow: operator.pow}
UNARY = {ast.UAdd: operator.pos, ast.USub: operator.neg}
FUNCTIONS = {"abs": abs, **{name: getattr(math, name)
                           for name in ("exp", "log", "sqrt", "sin", "cos", "tanh")}}
MODES = {"set": lambda _, value: value, "add": operator.add, "scale": operator.mul}


def number(value, context):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{context} must be a finite real number")
    return float(value)


def table(value, context, *, fields=None):
    if not isinstance(value, dict) or not value:
        raise ValueError(f"{context} must be a nonempty table")
    if fields is not None and (unknown := value.keys() - fields):
        raise ValueError(f"{context}: unknown fields {', '.join(sorted(unknown))}")
    return value


def expression(source, symbols):
    if not isinstance(source, str) or len(source) > 1024:
        raise ValueError("expressions must be text of at most 1024 characters")
    try:
        parsed = ast.parse(source.strip(), mode="eval")
    except SyntaxError as error:
        raise ValueError(f"invalid expression: {source!r}") from error
    if sum(1 for _ in ast.walk(parsed)) > 120:
        raise ValueError("expression is too complex")

    def read(node):
        if isinstance(node, ast.Constant):
            constant = number(node.value, "expression constant")
            return lambda _: constant
        if isinstance(node, ast.Name) and node.id in symbols:
            return lambda values: values[node.id]
        if isinstance(node, ast.BinOp) and type(node.op) in BINARY:
            fn = BINARY[type(node.op)]
            left, right = read(node.left), read(node.right)
            return lambda values: fn(left(values), right(values))
        if isinstance(node, ast.UnaryOp) and type(node.op) in UNARY:
            fn, item = UNARY[type(node.op)], read(node.operand)
            return lambda values: fn(item(values))
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id in FUNCTIONS and len(node.args) == 1 and not node.keywords):
            fn, item = FUNCTIONS[node.func.id], read(node.args[0])
            return lambda values: fn(item(values))
        raise ValueError(f"unsupported expression or unknown symbol: {ast.unparse(node)}")

    function = read(parsed.body)
    return lambda values: number(function(values), f"expression {source!r}")


@dataclass(frozen=True)
class Modulation:
    parameter: str
    start_s: float
    stop_s: float
    operation: str
    value: float


def read_modulations(items, parameters, end):
    if not isinstance(items, list) or len(items) > 1000:
        raise ValueError("modulations must be a list with at most 1000 entries")
    mods = []
    for item in items:
        fields = set(Modulation.__dataclass_fields__)
        table(item, "modulation", fields=fields)
        if item.keys() != fields:
            raise ValueError("a modulation requires parameter, start_s, stop_s, operation and value")
        parameter, op = item["parameter"], item["operation"]
        if not isinstance(parameter, str) or parameter not in parameters:
            raise ValueError("modulation parameter must name a declared parameter")
        if not isinstance(op, str) or op not in MODES:
            raise ValueError("modulation operation must be set, add or scale")
        start = number(item["start_s"], "modulation start_s")
        stop = number(item["stop_s"], "modulation stop_s")
        if not 0 <= start < stop <= end:
            raise ValueError("modulation must satisfy 0 <= start_s < stop_s <= time.stop_s")
        mods.append(Modulation(parameter, start, stop, op, number(item["value"], "modulation value")))
    mods.sort(key=lambda item: (item.parameter, item.start_s))
    for left, right in zip(mods, mods[1:]):
        if left.parameter == right.parameter and right.start_s < left.stop_s:
            raise ValueError(f"overlapping modulations for {left.parameter}")
    return mods


@dataclass(frozen=True)
class Reporter:
    source: object
    tau: float
    index: int | None


@dataclass(frozen=True)
class Detector:
    unit: str
    operator: str
    noise_sd: float | None
    weights: dict
    delay: Decimal
    exposure: Decimal
    index: int | None


class Acquisition:
    """Source response, detector mixing and a causal shutter."""

    def __init__(self, config, symbols, y):
        self.reporters = {}
        self.detectors = {}
        sources = table(config["reporters"], "reporters")
        channels = table(config["observations"], "observations")
        if len(sources) > 128 or len(channels) > 128:
            raise ValueError("ode-v2 supports at most 128 reporters and 128 detectors")
        for name, spec in sources.items():
            table(spec, name, fields={"expression", "response_s", "initial"})
            tau = number(spec.get("response_s", 0), f"{name} response_s")
            if tau < 0 or (tau > 0) != ("initial" in spec):
                raise ValueError(f"{name}: a positive response_s requires an explicit initial value")
            self.reporters[name] = Reporter(expression(spec["expression"], symbols), tau, len(y) if tau else None)
            if tau:
                y.append(number(spec["initial"], f"{name} initial"))
        self.exposures = {}
        for name, spec in channels.items():
            table(spec, name, fields={"unit", "operator", "noise_sd", "mix", "delay_s", "exposure_s"})
            unit, kind, sd = channel_metadata(name, spec)
            weights = {k: number(v, f"{name}/{k} mix coefficient")
                       for k, v in table(spec["mix"], f"{name} mix", fields=sources.keys()).items()}
            if not any(weights.values()):
                raise ValueError(f"{name}: at least one mix coefficient must be nonzero")
            delay = number(spec.get("delay_s", 0), f"{name} delay_s")
            exposure = number(spec.get("exposure_s", 0), f"{name} exposure_s")
            if min(delay, exposure) < 0:
                raise ValueError(f"{name}: delay_s and exposure_s must be nonnegative")
            detector = Detector(unit, kind, sd, weights, Decimal(str(delay)), Decimal(str(exposure)),
                                len(y) if exposure else None)
            self.detectors[name] = detector
            if exposure:
                self.exposures[name] = detector
                y.append(0.0)

    def response(self, name, context, values):
        reporter = self.reporters[name]
        return reporter.source(context) if reporter.index is None else values[reporter.index]

    def signal(self, detector, context, values):
        return math.fsum(weight * self.response(name, context, values)
                         for name, weight in detector.weights.items() if weight)

    def derivatives(self, context, values):
        response = [(reporter.source(context) - values[reporter.index]) / reporter.tau
                    for reporter in self.reporters.values() if reporter.tau]
        return response + [self.signal(d, context, values) for d in self.exposures.values()]

    def observe(self, times, edges, solutions, at_time):
        def segment(time):
            return min(bisect_right(edges, time) - 1, len(solutions) - 1)

        def average(detector, left, right, exposure):
            if left >= right:
                raise ValueError("exposure is below floating-point clock resolution")
            first, last = segment(left), segment(right)
            areas = []
            for i in range(first, last + 1):
                l, r = max(left, edges[i]), min(right, edges[i + 1])
                if l == r:
                    continue
                solve = solutions[i]
                before, after = solve(l)[detector.index], solve(r)[detector.index]
                area = after - before
                if abs(area) <= 1e-7 * max(abs(before), abs(after)):
                    # Avoid cancellation when shutter times nearly coincide on a long trajectory.
                    def integrand(t):
                        values = solve(t)
                        return self.signal(detector, at_time(t, values), values)
                    area, _ = quad(integrand, l, r, epsabs=ATOL, epsrel=RTOL)
                areas.append(area)
            return math.fsum(areas) / exposure

        observed = {}
        for name, detector in self.detectors.items():
            complete_times, readings = [], []
            for time in times:
                end = Decimal(str(time)) - detector.delay
                start = end - detector.exposure
                if start < 0:
                    continue
                endpoint = float(end)
                if detector.index is None:
                    last = solutions[segment(endpoint)](endpoint)
                    value = self.signal(detector, at_time(endpoint, last), last)
                else:
                    value = average(detector, float(start), endpoint, float(detector.exposure))
                readings.append(value)
                complete_times.append(time)
            if len(complete_times) <= int(detector.operator == "interval_rate"):
                raise ValueError(f"{name}: acquisition leaves no samples on the declared clock")
            observed[name] = observation(name, detector.unit, detector.operator, detector.noise_sd,
                                         np.asarray(complete_times), readings)
        return observed


def channel_metadata(name, channel):
    unit = channel["unit"]
    kind = channel.get("operator", "point")
    if not isinstance(unit, str) or not unit.strip() or kind not in ("point", "interval_rate"):
        raise ValueError(f"{name} requires a unit and a point or interval_rate operator")
    sd = number(channel["noise_sd"], f"{name} noise_sd") if "noise_sd" in channel else None
    if sd is not None and sd <= 0:
        raise ValueError(f"{name} noise_sd must be positive")
    return unit, kind, sd


def observation(name, unit, kind, sd, times, values):
    rate = kind == "interval_rate"
    series = (np.diff(values) / np.diff(times)).tolist() if rate else list(map(float, values))
    for value in series:
        number(value, f"{name} observation")
    return {"unit": unit, "operator": kind, "noise_sd": sd,
            "time_s": times[int(rate):].tolist(), "values": series}


def simulate(config):
    """Solve a portable TOML system and return observed channels in declared units."""
    table(config, "configuration")
    schema = config.get("schema")
    fields = {"schema", "time", "states", "parameters", "equations", "observations", "modulations"}
    if schema == ACQUISITION_SCHEMA:
        fields.add("reporters")
    table(config, "configuration", fields=fields)
    if schema not in SCHEMAS:
        raise ValueError(f"expected schema {SCHEMA} or {ACQUISITION_SCHEMA}")
    clock = table(config["time"], "time", fields={"stop_s", "sample_s"})
    end = number(clock["stop_s"], "stop_s")
    step = number(clock["sample_s"], "sample_s")
    if not 0 < step <= end or end / step > 100000:
        raise ValueError("require 0 < sample_s <= stop_s and at most 100001 samples")
    # Form decimal grid multiples before casting so a declared 0.9 boundary
    # matches the third step of a 0.3-second clock.
    decimal_step, decimal_end = Decimal(str(step)), Decimal(str(end))
    samples = int(decimal_end // decimal_step) + 1
    times = np.fromiter((float(index * decimal_step) for index in range(samples)), dtype=float)
    times = np.append(times[times < end], end)
    states = table(config["states"], "states")
    parameters = config.get("parameters", {})
    if not isinstance(parameters, dict) or len(states) + len(parameters) > 128:
        raise ValueError("parameters must be a table and states plus parameters must be <= 128")
    if states.keys() & parameters.keys():
        raise ValueError("states and parameter names must be distinct")
    symbols = states.keys() | parameters.keys()
    if any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) or keyword.iskeyword(name)
           for name in symbols):
        raise ValueError("state and parameter names must be ASCII identifiers")
    if symbols & (FUNCTIONS.keys() | {"t"}):
        raise ValueError("t and expression function names are reserved")
    y = [number(value, name) for name, value in states.items()]
    parameters = {name: number(value, name) for name, value in parameters.items()}
    eq = table(config["equations"], "equations", fields=states.keys())
    if eq.keys() != states.keys():
        raise ValueError("each state requires exactly one derivative equation")
    symbols = symbols | {"t"}
    derivatives = [expression(eq[name], symbols) for name in states]
    modulations = read_modulations(config.get("modulations", []), parameters, end)
    channels = {}
    acquisition = Acquisition(config, symbols, y) if schema == ACQUISITION_SCHEMA else None
    if acquisition is None:
        for name, channel in table(config["observations"], "observations").items():
            table(channel, name, fields={"expression", "unit", "operator", "noise_sd"})
            channels[name] = (*channel_metadata(name, channel), expression(channel["expression"], symbols))

    def parameters_at(time):
        active = parameters.copy()
        for item in modulations:
            if item.start_s <= time < item.stop_s:
                active[item.parameter] = number(MODES[item.operation](parameters[item.parameter], item.value),
                                                f"modulated {item.parameter}")
        return active

    def context(time, values, active):
        return {**active, "t": float(time), **dict(zip(states, map(float, values)))}

    edges = sorted({0.0, end, *(m.start_s for m in modulations), *(m.stop_s for m in modulations)})
    traces = np.empty((len(times), len(states)))
    solutions = []
    for left, right in zip(edges, edges[1:]):
        active = parameters_at(left)

        def rhs(time, values):
            current = context(time, values, active)
            rates = [equation(current) for equation in derivatives]
            if acquisition is not None:
                rates.extend(acquisition.derivatives(current, values))
            return rates

        solution = solve_ivp(rhs, (left, right), y, dense_output=True,
                             method="DOP853", rtol=RTOL, atol=ATOL)
        if not solution.success:
            raise ValueError(f"ODE integration failed: {solution.message}")
        selected = (times >= left) & ((times < right) if right < end else (times <= right))
        if np.any(selected):
            traces[selected] = solution.sol(times[selected])[:len(states)].T
        if acquisition is not None:
            solutions.append(solution.sol)
        y = solution.y[:, -1]
        if acquisition is not None:
            y = y.copy()
            for detector in acquisition.exposures.values():
                y[detector.index] = 0.0
    if not np.isfinite(traces).all():
        raise ValueError("ODE trajectory contains non-finite values")

    if acquisition is not None:
        observed = acquisition.observe(times, edges, solutions,
                                       lambda t, y: context(t, y, parameters_at(t)))
    else:
        values = {name: [] for name in channels}
        for time, state in zip(times, traces):
            current = context(time, state, parameters_at(float(time)))
            for name, (_, _, _, measure) in channels.items():
                values[name].append(measure(current))
        observed = {name: observation(name, unit, kind, sd, times, values[name])
                    for name, (unit, kind, sd, _) in channels.items()}
    return {"schema": schema, "time_s": times.tolist(), "states": dict(zip(states, traces.T.tolist())),
            "observations": observed,
            "solver": {"method": "DOP853", "rtol": RTOL, "atol": ATOL, "scipy": scipy_version}}


def compare(left, right):
    """Compare aligned expected observations without assigning mechanism probabilities."""
    a, b = simulate(left)["observations"], simulate(right)["observations"]
    if a.keys() != b.keys():
        raise ValueError("portable comparison requires the same observed channel names")
    results = {}
    for name, x in a.items():
        y = b[name]
        if any(x[key] != y[key] for key in ("unit", "operator", "time_s", "noise_sd")):
            raise ValueError(f"{name}: observation units, operators, times and noise scales must match")
        residual = np.asarray(x["values"]) - np.asarray(y["values"])
        if not len(residual):
            raise ValueError(f"{name}: comparison requires at least one sample")
        rmse = math.hypot(*residual) / math.sqrt(len(residual))
        sd = x["noise_sd"]
        results[name] = {"unit": x["unit"], "samples": len(residual), "rmse": number(rmse, "RMSE"),
                         "max_absolute_difference": number(float(np.max(abs(residual))), "maximum deviation"),
                         "rmse_in_noise_sd": number(rmse / sd, "scaled RMSE") if sd is not None else None}
    schema = ACQUISITION_SCHEMA if ACQUISITION_SCHEMA in (left["schema"], right["schema"]) else SCHEMA
    return {"schema": schema, "comparison": "deterministic expected observations", "channels": results}
