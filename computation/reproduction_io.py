"""Output paths and public input identity for the numerical reproduction."""
import hashlib
import os
from pathlib import Path

AMIRI_WORKBOOK_NAME = "41467_2023_43612_MOESM8_ESM.xlsx"
AMIRI_WORKBOOK_URL = (
    "https://www.ebi.ac.uk/biostudies/files/S-EPMC10695988/" + AMIRI_WORKBOOK_NAME
)
AMIRI_SHA256 = "dcf21574fbcd5318d1581b3e1c80bd77893c9214dd871a71f174bda49f23b907"


def output_directory(script: Path) -> Path:
    relative_parent = script.resolve().parent.relative_to(Path(__file__).resolve().parent)
    output = Path(os.environ["PROTRUSION_REPRO_OUTPUT"]) / relative_parent
    output.mkdir(parents=True, exist_ok=True)
    return output


def verify_workbook(contents: bytes) -> None:
    digest = hashlib.sha256(contents).hexdigest()
    if digest != AMIRI_SHA256:
        raise ValueError(f"Unexpected Amiri workbook SHA-256: {digest}")


def amiri_workbook() -> Path:
    workbook = Path(os.environ["PROTRUSION_AMIRI_WORKBOOK"])
    verify_workbook(workbook.read_bytes())
    return workbook
