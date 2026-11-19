"""Fetch and verify one source dataset and export plate-level observation rows."""

import argparse
from hashlib import sha256
from pathlib import Path
import subprocess
import tempfile
from urllib.request import urlopen

DOI = "10.1371/journal.pcbi.1014472"
COMMIT = "6161c370d7a616dde35a12f2a83fff4ee6e713aa"
SOURCE = f"https://raw.githubusercontent.com/snaketron/cellmig_paper_data/{COMMIT}/data"
INPUTS = {
    "calibration": ("dataset_1.RData", "489f8918e117753712d520b03ec3a855e96a76c262eaad07fb499b4e008b9d55"),
    "heldout": ("dataset_2.RData", "04323ec8db163b87bb2904a34a808a651a70193b488a5832812917f4d0c7c9bc"),
}


def extract(split, original):
    name, expected = INPUTS[split]
    if original is None:
        with urlopen(f"{SOURCE}/{name}", timeout=45) as source:
            raw = source.read(2_000_000)
    else:
        raw = original.read_bytes()
    if sha256(raw).hexdigest() != expected:
        raise ValueError(f"source SHA-256 mismatch for {name}")
    with tempfile.TemporaryDirectory() as temporary:
        here = Path(temporary) / name
        here.write_bytes(raw)
        return subprocess.check_output(
            ["Rscript", "--vanilla", str(Path(__file__).with_name("extract.R").resolve()), str(here), split],
            timeout=45,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", choices=tuple(INPUTS), required=True)
    parser.add_argument("--source", type=Path, help="pinned original RData file (optional)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = extract(args.split, args.source)
        with args.output.open("xb") as stream:
            stream.write(result)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
