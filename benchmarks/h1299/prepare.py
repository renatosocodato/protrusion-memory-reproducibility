"""Export paired H1299 observations from the pinned PLOS supplementary XLSB."""

import argparse
import csv
from datetime import date as calendar_date
from hashlib import sha256
from pathlib import Path
import tempfile
import tomllib
from urllib.request import urlopen

from pyxlsb import open_workbook

ARTICLE = "10.1371/journal.pone.0135204"
SOURCE = f"https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/{ARTICLE}.s001"
SOURCE_SHA256 = "5d6cda93f76ca04c02f71459f4fab87365b99d2b6b9f8a420ecc6fbeb1d7d4d3"
AREA_PER_PIXEL_UM2 = 0.21 ** 2
HEADERS = ("preparation", "unit", "frame", "actin", "protrusion")


def download(source):
    if source is not None:
        return source.read_bytes()
    with urlopen(SOURCE, timeout=60) as response:
        return response.read(12_000_000)


def extract(original, dates):
    if sha256(original).hexdigest() != SOURCE_SHA256:
        raise ValueError("public XLSB SHA-256 mismatch")
    with tempfile.TemporaryDirectory() as temp:
        workbook = Path(temp) / "h1299.xlsb"
        workbook.write_bytes(original)
        with open_workbook(workbook) as source, source.get_sheet("Features") as features, \
                source.get_sheet("Processes") as motion, source.get_sheet("Index") as indices:
            streams = [iter(sheet.rows()) for sheet in (features, motion, indices)]
            fhead, mhead, ihead = [[c.v for c in next(stream)] for stream in streams]
            actin = fhead.index("Cell Mean LifeAct")
            protrusion = mhead.index("Protrusion Pixels")
            date, cell, frame = (ihead.index(key) for key in ("Date", "Cell Trace", "Frame Number"))
            exported, seen = [], set()
            for f, m, i in zip(*streams, strict=True):
                if not f[0].r == m[0].r == i[0].r:
                    raise ValueError("source worksheet rows do not align")
                source_date = i[date].v
                preparation = str(int(source_date))
                if source_date != int(source_date):
                    raise ValueError("experimental date is not an integer YYYYMMDD label")
                calendar_date.fromisoformat(preparation)
                key = (preparation, str(i[cell].v), int(i[frame].v))
                if key in seen or i[frame].v != key[2]:
                    raise ValueError("invalid or duplicate source frame")
                seen.add(key)
                if preparation not in dates:
                    continue
                target = m[protrusion].v
                exported.append((*key, f[actin].v,
                                 target * AREA_PER_PIXEL_UM2 if target is not None else None))
    if (len(seen), len({k[0] for k in seen}), len({k[:2] for k in seen})) != (6419, 19, 122):
        raise ValueError("source observation, preparation or track counts do not match the study")
    return exported


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", choices=("calibration", "heldout"), required=True)
    parser.add_argument("--source", type=Path, help="original public XLSB file")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with Path(__file__).with_name("forecast.toml").open("rb") as stream:
        dates = set(tomllib.load(stream)["split"][args.split])
    data = extract(download(args.source), dates)
    with args.output.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(HEADERS)
        writer.writerows(data)


if __name__ == "__main__":
    main()
