"""Test selective mechanism decisions on independent synthetic contrast families."""

from dataclasses import dataclass
import json

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

MIN_EFFECT = .08
ABSTAIN = -1


@dataclass
class Panel:
    h: np.ndarray
    c: np.ndarray
    coarse: np.ndarray
    width: np.ndarray
    truth: np.ndarray
    twin: np.ndarray
    qualified: np.ndarray

    def prior_decision(self):
        return decide(self.h, self.coarse + self.width)

    def control_decision(self, require_qualified):
        rule = decide(self.h - self.c, self.width + .018)
        return np.where(self.qualified | (not require_qualified), rule, ABSTAIN)

    def features(self):
        # The independent measurement and the true pathway are unavailable at acquisition.
        return np.column_stack((self.h, self.coarse, self.width,
                                np.abs(self.h - MIN_EFFECT) - self.coarse - self.width))


def decide(effect, width):
    return np.select([effect - width > MIN_EFFECT, effect + width <= MIN_EFFECT], [1, 0], default=ABSTAIN)


def panel(seed, family_count, domain):
    rng = np.random.default_rng(seed)
    # Independent family contrasts. The last family class has two exactly matched observed twins.
    easy = rng.integers(0, 2, family_count)
    a = np.concatenate((easy * rng.uniform(.20, .34, family_count), np.zeros(family_count),
                        rng.uniform(.21, .34, family_count)))
    width = rng.uniform(.010, .024 if domain == "development" else .038, family_count)
    twin_total = a[2 * family_count:] + rng.uniform(.005, .025, family_count)
    true_total = np.concatenate((a[:family_count] + rng.uniform(.006, .020, family_count),
                                 twin_total, twin_total))
    observation_noise = rng.uniform(-width / 2, width / 2, (2, family_count))
    h = true_total + np.concatenate((observation_noise[0], observation_noise[1], observation_noise[1]))
    coarse = np.repeat([.035, .40, .40], family_count)
    width = np.tile(width, 3)
    truth = (a > MIN_EFFECT).astype(int)
    twins = coarse > .1
    c = true_total - a + rng.uniform(-.009, .009, 3 * family_count)
    qualified = np.ones_like(twins)
    if domain == "nontransfer":
        # The actin-disabled preparation omits the extra direct response in the intact direct twin.
        c[family_count:2 * family_count] -= a[2 * family_count:]
        qualified[family_count:] = False  # a separate mechanical witness rules out transport in this stratum
    out = Panel(h, c, coarse, width, truth, twins, qualified)
    np.testing.assert_array_equal(out.features()[family_count:2 * family_count],
                                  out.features()[2 * family_count:])
    return out


def fit_logistic(x, y, weight):
    standard = np.maximum(x.std(axis=0), .001)
    center = x.mean(axis=0)
    design = np.column_stack((np.ones(len(x)), (x - center) / standard))

    def loss(beta):
        z = design @ beta
        residual = weight * (expit(z) - y)
        return (np.mean(weight * (np.logaddexp(0, z) - y * z)) + .001 * (beta @ beta),
                (design.T @ residual) / len(y) + .002 * beta)

    optimized = minimize(loss, np.zeros(design.shape[1]), jac=True, method="L-BFGS-B")
    if not optimized.success:
        raise RuntimeError(optimized.message)
    return lambda x: expit(optimized.x[0] + ((x - center) / standard) @ optimized.x[1:])


def summarize(decision, requested, panel):
    decided = decision != ABSTAIN
    wrong = decided & (decision != panel.truth)
    direct_twins = panel.twin & (panel.truth == 0)
    mediated_twins = panel.twin & (panel.truth == 1)
    return {
        "measurements_per_case": float(np.mean(requested)),
        "coverage": float(decided.mean()),
        "error_among_calls": float(wrong.sum() / decided.sum()) if decided.any() else None,
        "false_mediation_among_direct_twins": float(np.mean(decision[direct_twins] == 1)),
        "coverage_direct_twins": float(decided[direct_twins].mean()),
        "coverage_mediated_twins": float(decided[mediated_twins].mean()),
    }


def benchmark():
    train = panel(72913, 1200, "development")
    acquisition = fit_logistic(train.features(), (train.prior_decision() == ABSTAIN)
                               & (train.control_decision(True) == train.truth), np.ones(len(train.h)))
    # Reweight the exactly identical ambiguous observations to 95% mediation in development.
    # Held-out families are balanced. A prevalence shortcut cannot recover the hidden pathway.
    prior_weight = np.where(train.twin, np.where(train.truth == 1, 1.9, .1), 1)
    passive = fit_logistic(train.features(), train.truth, prior_weight)
    result = {"data_kind": "synthetic contrast-level surrogate", "cost_per_mechanical_control": 1,
              "development_families": 2400, "test_families_per_domain": 2400,
              "minimum_effect_um": MIN_EFFECT, "domains": {}}
    for domain, seed in (("transfer", 72914), ("nontransfer", 72915)):
        test = panel(seed, 1200, domain)
        prior = test.prior_decision()
        scores = passive(test.features())
        passive_call = np.where(scores >= .75, 1, np.where(scores <= .25, 0, ABSTAIN))
        analytic_request = prior == ABSTAIN
        learned_request = (acquisition(test.features()) >= .5) & analytic_request
        requests = {"analytic": analytic_request, "learned": learned_request,
                    "always_calibrate": np.ones_like(test.twin)}
        # Use the mechanical witness before purchasing the replicate/control in every qualified policy.
        out = {"passive_prevalence": summarize(passive_call, 0, test),
               "analytic_no_witness": summarize(np.where(analytic_request, test.control_decision(False), prior),
                                                analytic_request, test)}
        for name, request in requests.items():
            request = request & test.qualified
            if name == "always_calibrate":
                decision = test.control_decision(True)
            else:
                decision = np.where(request, test.control_decision(True), prior)
            out[name] = summarize(decision, request, test)
        result["domains"][domain] = out
    return result


if __name__ == "__main__":
    print(json.dumps(benchmark(), indent=2, allow_nan=False))
