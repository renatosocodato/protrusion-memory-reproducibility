"""Package a committed source tree and its regenerated public paper files."""

import argparse
from hashlib import sha256
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
DOCUMENTS = ("manuscript.pdf", "manuscript.docx",
             "supporting_information.pdf", "supporting_information.docx")


def bundle(rendered, output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise ValueError("package the paper from a clean, committed source checkout")
    with (ROOT / "pyproject.toml").open("rb") as project:
        version = tomllib.load(project)["project"]["version"]
    name = f"protrusion-memory-paper-{version}"
    output.mkdir(parents=True, exist_ok=True)
    if (output / f"{name}.zip").exists():
        raise FileExistsError(f"bundle already exists: {name}.zip")
    with tempfile.TemporaryDirectory() as temporary:
        workspace = Path(temporary)
        archive = subprocess.check_output(["git", "archive", f"--prefix={name}/", "HEAD"], cwd=ROOT)
        with tarfile.open(fileobj=io.BytesIO(archive)) as source:
            source.extractall(workspace, filter="data")
        root = workspace / name
        products = root / "rendered"
        products.mkdir()
        for filename in DOCUMENTS:
            shutil.copyfile(rendered / filename, products / filename)
        figures = products / "figures"
        figures.mkdir()
        for source_file in sorted((rendered / "figures").rglob("*")):
            if source_file.suffix not in {".pdf", ".svg", ".png"}:
                continue
            path = figures / source_file.relative_to(rendered / "figures")
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_file, path)
        manifest = {path.relative_to(root).as_posix(): sha256(path.read_bytes()).hexdigest()
                    for path in sorted(root.rglob("*")) if path.is_file()}
        (root / "paper-bundle-sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")
        return shutil.make_archive(str(output / name), "zip", root_dir=workspace, base_dir=name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rendered", type=Path, required=True,
                        help="validated article/submission_biophysj output directory")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(bundle(args.rendered.resolve(), args.output.resolve()))


if __name__ == "__main__":
    main()
