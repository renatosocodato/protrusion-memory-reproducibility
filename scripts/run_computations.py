#!/usr/bin/env python3
"""Rerun the synthetic calculations and the published source-data summaries."""
import argparse
from collections import defaultdict
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from urllib.request import urlopen

COMPUTATION = Path(__file__).resolve().parents[1] / "computation"
sys.dont_write_bytecode = True
sys.path.insert(0, str(COMPUTATION))
from reproduction_io import AMIRI_WORKBOOK_NAME, AMIRI_WORKBOOK_URL, verify_workbook

MAIN = {
    "joint_fret_actin_mechanics": (
        "run_joint_observation",
        "verify_twin_symbolic",
        "independent_gate_readout",
    ),
    "frontier_conformation_chaos": ("run_frontier", "feedback_boundary"),
    "autonomous_ensemble_chaos": ("run_autonomous", "run_dop853_check"),
    "wasl_omics_gate": ("run_wasl_sensitivity",),
}

SUPPLEMENT = {
    "mechanochemical_disorder": ("run_mechanochemical",),
    "autonomous_ensemble_chaos": (
        "plot_autonomous",
        "formalism_note/check_identities",
        "formalism_note/check_symmetry",
    ),
    "molecular_organization_protrusion": (
        "run_comparison",
        "paired_uncertainty",
        "exploratory_matched_recruitment_v2",
    ),
    "spatial_molecular_protrusion": (
        "run_spatial",
        "second_challenge",
        "third_challenge",
        "fourth_challenge",
        "check_nonlinear_formula",
        "render_joint_challenge",
        "render_harmonics",
    ),
    "conformational_ensemble": (
        "run_challenge",
        "replacement_holdout",
        "filter_twin",
        "render_validated_comparison",
    ),
}

EXPLORATORY = {"wasl_omics_gate": ("run_wasl_dynamics", "simulate_wasl_measurement")}

AMIRI = {
    "mechanochemical_disorder": ("inspect_amiri", "plot_comparison"),
    "wasl_paired_reconciliation": ("calibration_reference_screen",),
}

THREAD_LIMITS = dict.fromkeys(
    (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "NUMBA_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ),
    "1",
)


def workbook_path(output: Path, local: Path | None) -> Path:
    if local is not None:
        verify_workbook(local.read_bytes())
        return local.resolve()

    with urlopen(AMIRI_WORKBOOK_URL, timeout=30) as response:
        contents = response.read(1_000_000)
    verify_workbook(contents)
    cached = output / "inputs" / AMIRI_WORKBOOK_NAME
    cached.parent.mkdir()
    cached.write_bytes(contents)
    return cached


def run_steps(suites, output: Path, env: dict[str, str]) -> dict[str, float]:
    elapsed = {}
    for suite in suites:
        for experiment, scripts in suite.items():
            for name in scripts:
                relative = Path(experiment) / (name + ".py")
                log = output / "logs" / relative.with_suffix(".log")
                log.parent.mkdir(parents=True, exist_ok=True)
                print(f"Running {relative}", flush=True)
                start = time.monotonic()
                try:
                    with log.open("w") as stream:
                        subprocess.run(
                            [sys.executable, str(COMPUTATION / relative)],
                            cwd=output,
                            env=env,
                            stdout=stream,
                            stderr=subprocess.STDOUT,
                            check=True,
                            timeout=150,
                        )
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                    raise RuntimeError(f"{relative} failed. See {log}") from error
                elapsed[str(relative)] = round(time.monotonic() - start, 2)
    return elapsed


def same_value(observed, reference, tolerance):
    if isinstance(reference, list):
        return isinstance(observed, list) and len(observed) == len(reference) and all(
            same_value(actual, expected, tolerance)
            for actual, expected in zip(observed, reference)
        )
    if type(reference) is float:
        return math.isclose(observed, reference, rel_tol=0.0, abs_tol=tolerance)
    return type(observed) is type(reference) and observed == reference


def verify_results(output: Path, groups: list[str]) -> dict[str, int]:
    checks = json.loads((COMPUTATION / "reference_checks.json").read_text())
    documents = {}
    counts = defaultdict(int)
    for group in groups:
        for check in checks[group]:
            path = check["file"]
            if path not in documents:
                documents[path] = json.loads((output / path).read_text())
            for field, expected in check["values"].items():
                actual = documents[path]
                for key in field.split("/"):
                    actual = actual[int(key)] if isinstance(actual, list) else actual[key]
                if not same_value(actual, expected, check.get("absolute_tolerance", 1e-9)):
                    raise ValueError(f"{path} [{field}] produced {actual!r}; expected {expected!r}")
                counts[group] += 1
    return dict(counts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new or empty output directory")
    parser.add_argument("--scope", choices=("main", "supplement", "all"), default="main")
    public = parser.add_mutually_exclusive_group()
    public.add_argument("--amiri-workbook", type=Path, help="existing, unmodified public Source Data XLSX")
    public.add_argument("--skip-public", action="store_true", help="run synthetic and analytic calculations offline")
    parser.add_argument("--include-exploratory", action="store_true",
                        help="also run abundance scenarios retained outside the article")
    args = parser.parse_args()
    if args.scope == "main" and args.amiri_workbook:
        parser.error("--amiri-workbook requires --scope supplement or all")

    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("--output must be a new or empty directory")
    output.mkdir(parents=True, exist_ok=True)

    env = {
        **os.environ,
        **THREAD_LIMITS,
        "PYTHONPATH": str(COMPUTATION),
        "PYTHONDONTWRITEBYTECODE": "1",
        "MPLCONFIGDIR": str(output / ".matplotlib"),
        "PROTRUSION_REPRO_OUTPUT": str(output),
        "PROTRUSION_REPRO_EXPLORATORY": "1" if args.include_exploratory else "0",
    }
    suites = [MAIN]
    groups = ["main"]
    if args.scope != "main":
        suites.append(SUPPLEMENT)
        groups.append("supplement")
        if not args.skip_public:
            try:
                env["PROTRUSION_AMIRI_WORKBOOK"] = str(workbook_path(output, args.amiri_workbook))
            except (OSError, ValueError) as error:
                raise SystemExit(
                    f"Unable to verify the public Amiri workbook: {error}. "
                    "Use --amiri-workbook with the publisher's workbook or --skip-public for an offline run."
                ) from error
            suites.append(AMIRI)
            groups.append("amiri")

    if args.include_exploratory:
        suites.append(EXPLORATORY)
        groups.append("exploratory")

    elapsed = run_steps(suites, output, env)
    counts = verify_results(output, groups)
    report = {"scope": args.scope, "reference_checks_passed": counts, "script_seconds": elapsed}
    (output / "reproduction_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Passed {sum(counts.values())} reference checks. Results: {output}")


if __name__ == "__main__":
    main()
