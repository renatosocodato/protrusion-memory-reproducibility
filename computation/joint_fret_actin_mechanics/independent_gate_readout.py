"""Sensitivity of the published synthetic twin to an additional calibrated readout."""
import json
from pathlib import Path
from reproduction_io import output_directory

import numpy as np
from scipy.linalg import expm
from scipy.special import ndtr
from scipy.integrate import solve_ivp

TAUS = (0, 10, 30)
CADENCES = (5, 15, 60)
MIXING_BOUNDS = (0.20, 0.50, 0.68)
NOISE = 0.03
SEED = 202611180301
DRAWS = 10000
CRITICAL = 1.959963984540054


def filtered_contrast(tau, order):
    """Propagate q_A, a and, if needed, the filtered q_A minus a."""
    dynamics = np.zeros((5, 5))
    dynamics[0, [0, 3]] = [-1 / 8, 1 / 8]
    dynamics[1, [0, 1, 4]] = [1 / 25, -1 / 25, 1 / 25]
    if tau:
        dynamics[2, :3] = np.array([1, -1, -1]) / tau
    step = expm(dynamics)
    state = np.zeros(3)
    contrast = [0.0]
    windows = [(20, 80), (140, 200)]
    if order == "vu":
        windows.reverse()
    for t in range(300):
        u = 0.3 * (windows[0][0] <= t < windows[0][1])
        v = 0.2 * (windows[1][0] <= t < windows[1][1])
        state = (step @ np.r_[state, u, v])[:3]
        contrast.append(state[2] if tau else state[0] - state[1])
    return np.array(contrast)


def independent_check():
    """Check filter dynamics by adaptive integration across pulse boundaries."""
    states = np.zeros(3)
    samples = [0.0]
    boundaries = (0, 20, 80, 140, 200, 300)
    for start, end in zip(boundaries, boundaries[1:]):
        u = 0.3 if start == 20 else 0.0
        v = 0.2 if start == 140 else 0.0

        def rhs(_, state):
            gate, actin, sensor = state
            return ((u - gate) / 8, (gate + v - actin) / 25,
                    (gate - actin - sensor) / 10)

        solution = solve_ivp(rhs, (start, end), states, method="DOP853",
                             t_eval=np.arange(start + 1, end + 1),
                             rtol=1e-11, atol=1e-13)
        if not solution.success:
            raise RuntimeError(solution.message)
        samples.extend(solution.y[2])
        states = solution.y[:, -1]
    return float(np.max(np.abs(samples - filtered_contrast(10, "uv"))))


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for tau in TAUS:
        contrasts = {order: filtered_contrast(tau, order) for order in ("uv", "vu")}
        for cadence in CADENCES:
            norms = [(np.linalg.norm(g[phase::cadence]), order, phase)
                     for order, g in contrasts.items() for phase in range(cadence)]
            smallest, order, phase = min(norms)
            largest = max(norms)[0]
            for bound in MIXING_BOUNDS:
                # Under the projection, closest means are at m_A=b and m_B=0.
                gap = (0.68 - bound) * smallest / NOISE
                # The thresholds also bound the probability of a wrong-family
                # declaration at the nearest endpoint. Between them we abstain.
                threshold = max(gap / 2, CRITICAL)
                distance = gap - threshold
                # Each draw represents an independently generated noisy trace
                # reduced to its sufficient scalar, under the known templates.
                noise = rng.standard_normal((2, DRAWS))
                empirical = np.mean(noise < distance, axis=1)
                wrong = np.mean(noise > threshold, axis=1)
                rows.append(dict(tau=tau, cadence=cadence, bound=bound, order=order,
                                 phase=phase, min_norm=smallest, max_norm=largest,
                                 gap_sigma=gap, correct=float(ndtr(distance)),
                                 wrong=float(ndtr(-threshold)),
                                 mc_correct=empirical.tolist(), mc_wrong=wrong.tolist()))
    assert abs(max(abs(filtered_contrast(0, "uv"))) * .68 - .126) < .001
    assert all(row["wrong"] <= .025000000001 for row in rows)
    assert all(abs(row["correct"] - x) < .02 for row in rows for x in row["mc_correct"])
    integration_error = independent_check()
    assert integration_error < 1e-10
    # At b=0.68 the closest means and the midpoint coincide.
    null = rng.standard_normal(DRAWS)
    result = dict(data_kind="synthetic", causal_use="observation design under fixed rival models",
                  seed=SEED, draws_per_family_and_setting=DRAWS, sigma=NOISE,
                  equal_mean_null_declared=float(np.mean(abs(null) > CRITICAL)),
                  equal_mean_null_theory=.05,
                  integration_max_abs_difference=integration_error, rows=rows)
    output = output_directory(Path(__file__)) / "independent_gate_readout.json"
    output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
