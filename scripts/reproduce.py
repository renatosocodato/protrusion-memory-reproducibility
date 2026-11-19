"""Rebuild the artwork and reading copies, optionally rerunning the analyses."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PAPER_FILES = (
    "AUTHOR_FRONT_MATTER.md", "SIGNIFICANCE.md", "bj_reading.latex",
    "build_submission.py",
)


def run(*command, env=None):
    subprocess.run([str(part) for part in command], cwd=ROOT, env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build")
    parser.add_argument("--artwork-only", action="store_true",
                        help="Render from supplied figure inputs without rerunning analyses.")
    parser.add_argument("--amiri-workbook", type=Path,
                        help="Use a local copy of the checksum-verified public workbook.")
    args = parser.parse_args()
    if args.artwork_only and args.amiri_workbook is not None:
        parser.error("--amiri-workbook is only used for full reproduction")
    build = args.output.resolve()
    if build.exists() and (not build.is_dir() or any(build.iterdir())):
        parser.error("--output must be a new or empty directory")
    required = ("Rscript", "pdfinfo", "pandoc", "xelatex")
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        parser.error(f"missing system dependencies: {', '.join(missing)}")
    manuscript = build / "article"
    submission = manuscript / "submission_biophysj"
    figures = submission / "figures"
    submission.mkdir(parents=True, exist_ok=True)

    environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                       VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1",
                       REPRO_DATA_DIR=str(ROOT / "figures" / "data"))
    if not args.artwork_only:
        output = build / "computations"
        command = [sys.executable, ROOT / "scripts/run_computations.py",
                   "--scope", "all", "--output", output]
        if args.amiri_workbook:
            command.extend(["--amiri-workbook", args.amiri_workbook.resolve()])
        run(*command, env=environment)
        data = build / "figure_inputs"
        run(sys.executable, ROOT / "scripts/regenerate_figure_inputs.py", output, data, env=environment)
        environment["REPRO_DATA_DIR"] = str(data)

    for name in ("render_main.R", "render_supplementary.R"):
        run("Rscript", "--vanilla", Path("figures") / name, figures, env=environment)
    run(sys.executable, ROOT / "scripts/check_figures.py", figures)
    for name in PAPER_FILES:
        shutil.copyfile(ROOT / "paper" / name, submission / name)
    for name in ("MANUSCRIPT_WORKING_DRAFT_EN.md", "check_callout_order.py"):
        shutil.copyfile(ROOT / "paper" / name, manuscript / name)
    run(sys.executable, submission / "build_submission.py", env=environment)
    print(f"PDF and Word reading copies: {submission}")


if __name__ == "__main__":
    main()
