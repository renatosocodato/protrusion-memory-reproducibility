"""Command-line interface for modeling and paper reproduction."""

import argparse
from dataclasses import asdict
from importlib.metadata import version
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tomllib

from .model import Model, Pulse, Reporter, compare_reporter, simulate_pair


def read_config(path):
    if path is None:
        return Model(), Reporter()
    with path.open("rb") as stream:
        sections = tomllib.load(stream)
    unknown = sections.keys() - {"model", "gate_pulse", "actin_pulse", "reporter"}
    if unknown:
        raise ValueError(f"unknown TOML section: {', '.join(sorted(unknown))}")
    defaults = Model()
    return Model(
        **sections.get("model", {}),
        gate=Pulse(**{**asdict(defaults.gate), **sections.get("gate_pulse", {})}),
        actin=Pulse(**{**asdict(defaults.actin), **sections.get("actin_pulse", {})}),
    ), Reporter(**sections.get("reporter", {}))


def simulation(model, order):
    a, b, optical = simulate_pair(model, order)
    return {
        "description": "Synthetic deviations in assigned model units",
        "order": order,
        "time_s": list(range(model.duration_s + 1)),
        "columns": ["gate", "actin", "deformation", "optical_precursor"],
        "A_fast_gate": [[*state, signal] for state, signal in zip(a.tolist(), optical.tolist())],
        "B_slow_gate": [[*state, state[0]] for state in b.tolist()],
    }


def reproduction(args, installed_version):
    root = args.bundle.resolve()
    with (root / "pyproject.toml").open("rb") as stream:
        project = tomllib.load(stream)["project"]
    if project["name"] != "protrusion-memory" or project["version"] != installed_version:
        raise ValueError(f"use a Protrusion Memory {installed_version} paper bundle or checkout")
    manifest = root / "paper-bundle-sha256.json"
    if not manifest.is_file() and not (root / ".git").exists():
        raise ValueError("paper bundle SHA-256 manifest is missing")
    if manifest.is_file():
        entries = json.loads(manifest.read_text(encoding="utf-8"))
        if "scripts/reproduce.py" not in entries:
            raise ValueError("paper bundle manifest omits the reconstruction script")
        for name, expected in entries.items():
            item = (root / name).resolve()
            if not item.is_relative_to(root) or hashlib.sha256(item.read_bytes()).hexdigest() != expected:
                raise ValueError(f"paper bundle checksum mismatch: {name}")
    runner = root / "scripts" / "reproduce.py"
    if not runner.is_file():
        raise ValueError("provide --bundle with the extracted paper bundle or full source checkout")
    command = [sys.executable, str(runner), "--output", str(args.output.resolve())]
    if args.artwork_only:
        command.append("--artwork-only")
    if args.amiri_workbook:
        command.extend(("--amiri-workbook", str(args.amiri_workbook.resolve())))
    return subprocess.run(command, cwd=root, check=False).returncode


def main(argv=None):
    installed_version = version("protrusion-memory")
    parser = argparse.ArgumentParser(prog="protrusion-memory", description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {installed_version}")
    commands = parser.add_subparsers(dest="command", required=True)
    sim = commands.add_parser("simulate", help="simulate the exact linear model pair")
    compare = commands.add_parser("compare", help="compare the calibrated reporter predictions")
    for command in (sim, compare):
        command.add_argument("--config", type=Path, help="optional TOML with model and reporter settings")
        command.add_argument("--output", type=Path, help="new JSON file (default: standard output)")
    sim.add_argument("--order", choices=("uv", "vu"), default="uv",
                     help="use configured windows (uv) or swap them (vu)")
    reproduce = commands.add_parser("reproduce", help="rebuild the frozen paper workflow")
    reproduce.add_argument("--target", choices=("paper",), default="paper")
    reproduce.add_argument("--bundle", type=Path, default=Path.cwd(),
                           help="extracted paper bundle or full checkout (default: current directory)")
    reproduce.add_argument("--output", type=Path, default=Path("build") / "paper")
    reproduce.add_argument("--artwork-only", action="store_true")
    reproduce.add_argument("--amiri-workbook", type=Path, help="original public workbook for offline use")
    args = parser.parse_args(argv)
    try:
        if args.command == "reproduce":
            return reproduction(args, installed_version)
        model, reporter = read_config(args.config)
        result = {"software_version": installed_version, "model": asdict(model)}
        if args.command == "simulate":
            result.update(simulation(model, args.order))
        else:
            result.update(compare_reporter(model, reporter))
            result["reporter"] = asdict(reporter)
            result["assumptions"] = ["known unit gain", "exact common actin subtraction",
                                     "independent Gaussian noise with fixed SD per sample",
                                     "single-record decision without multiplicity adjustment"]
        text = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.output is None:
            sys.stdout.write(text)
        else:
            with args.output.open("x", encoding="utf-8") as stream:
                stream.write(text)
    except (OSError, ValueError, TypeError, KeyError, ArithmeticError) as error:
        if isinstance(error, FileNotFoundError) and args.command == "reproduce":
            parser.error(f"{error}. Provide --bundle with the extracted paper bundle or full checkout.")
        parser.error(str(error))
    return 0
