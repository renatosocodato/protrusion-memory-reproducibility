#!/usr/bin/env python3
"""Check the leading nonlinear terms by direct periodic finite differences."""
import json
from pathlib import Path

import numpy as np

from reproduction_io import output_directory


def main():
    points = 64
    spacing = 32.0 / points
    mobility = 0.2
    beta = 1.0
    phase = 2 * np.pi * 2 * np.arange(points) / points + 0.3
    fundamental = np.cos(phase)
    third_harmonic = np.cos(3 * phase)
    cubed = fundamental**3
    lambda_2, lambda_6 = 4 * np.sin(np.pi * np.array([2, 6]) / points)**2 / spacing**2
    gamma = mobility * beta * lambda_6
    laplacian = (np.roll(cubed, -1) - 2 * cubed + np.roll(cubed, 1)) / spacing**2
    sources = {
        "D": mobility * beta * laplacian,
        "H": -gamma * (cubed - cubed.mean()),
    }
    expected = {"D": lambda_6 / (3 * lambda_2), "H": 1 / 3}
    observed = {}
    for model, source in sources.items():
        c1 = -2 * np.mean(source * fundamental)
        c3 = -2 * np.mean(source * third_harmonic)
        ratio = c3 / c1
        assert np.isclose(ratio, expected[model], rtol=0, atol=1e-12)
        observed[model] = {"C1": float(c1), "C3": float(c3), "Q": float(ratio)}
    result = {
        "method": "direct real-space central differences, independent of Fourier time step",
        "expected_D": float(expected["D"]),
        "expected_H": float(expected["H"]),
        "observed": observed,
    }
    output = output_directory(Path(__file__)) / "results_v1" / "nonlinear_formula_independent_check.json"
    output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
