#!/usr/bin/env python3
"""Check dimensions and formats of all regenerated publication figures."""

import argparse
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import xml.etree.ElementTree as ET
import zlib


DPI = 600
POINTS_PER_MM = 72 / 25.4
SVG_NAMESPACE = "{http://www.w3.org/2000/svg}"

# Physical dimensions in millimetres, in manuscript order.
PLATES = {
    "F01": (180, 162),
    "F02": (180, 187),
    "F03": (180, 167),
    "S01": (180, 66),
    "S02": (180, 63),
    "S03": (180, 117),
    "S04": (180, 151),
    "S05": (180, 131),
    "S06": (180, 121),
    "S07": (180, 74),
}
MAIN_PANELS = {
    "F01": [(180, 38)] + [(90, 63)] * 4,
    "F02": [(180, 38)] + [(90, 63)] * 4 + [(180, 37)],
    "F03": [(180, 38), (90, 63), (90, 63), (180, 63)],
}
SUPPLEMENTARY_GRIDS = {
    # Rows and relative column widths.
    "S01": (1, (1, 1, 1)),
    "S02": (1, (1, 1, 1)),
    "S03": (2, (1, 1, 1)),
    "S04": (2, (1, 1)),
    "S05": (2, (1, 1, 1.15)),
    "S06": (2, (1, 1)),
    "S07": (1, (1, 1)),
}


def expected_sizes():
    sizes = dict(PLATES)
    for figure, panels in MAIN_PANELS.items():
        for index, size in enumerate(panels):
            sizes[f"panels/{figure}{chr(ord('a') + index)}"] = size
    for figure, (rows, columns) in SUPPLEMENTARY_GRIDS.items():
        width, height = PLATES[figure]
        for index, column in enumerate(columns * rows):
            sizes[f"panels/{figure}{chr(ord('a') + index)}"] = (
                width * column / sum(columns), height / rows
            )
    return sizes


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_png(path, width, height):
    expected_pixels = tuple(int(mm * DPI / 25.4) for mm in (width, height))
    expected_density = round(DPI / 0.0254)
    with path.open("rb") as png:
        require(png.read(8) == b"\x89PNG\r\n\x1a\n", "invalid PNG signature")
        pixels = density = None
        image_data_seen = False
        while True:
            header = png.read(8)
            require(len(header) == 8, "truncated PNG")
            length, kind = struct.unpack(">I4s", header)
            payload = png.read(length)
            require(len(payload) == length, "truncated PNG chunk")
            crc = png.read(4)
            require(crc == struct.pack(">I", zlib.crc32(payload, zlib.crc32(kind))),
                    f"invalid {kind!r} checksum")
            if kind == b"IHDR":
                pixels = struct.unpack(">II", payload[:8])
            elif kind == b"pHYs":
                density = struct.unpack(">IIB", payload)
            elif kind == b"IDAT":
                image_data_seen = True
            elif kind == b"IEND":
                break
        require(image_data_seen, "missing PNG image data")
        require(pixels == expected_pixels, f"expected {expected_pixels} pixels, found {pixels}")
        require(density == (expected_density, expected_density, 1),
                f"expected {DPI} dpi in both directions")


def check_svg(path, width, height):
    svg = ET.parse(path).getroot()
    require(svg.tag == f"{SVG_NAMESPACE}svg", "invalid SVG root")
    viewbox = tuple(float(value) for value in svg.attrib["viewBox"].split())
    sizes = tuple(float(svg.attrib[attr].removesuffix("pt")) for attr in ("width", "height"))
    expected = (width * POINTS_PER_MM, height * POINTS_PER_MM)
    require(len(viewbox) == 4 and viewbox[:2] == (0, 0), "invalid SVG viewBox")
    for actual in (viewbox[2:], sizes):
        require(all(abs(value - target) < 0.02 for value, target in zip(actual, expected)),
                f"expected SVG size {expected} pt, found {actual}")
    tags = {element.tag.removeprefix(SVG_NAMESPACE) for element in svg.iter()}
    require("text" in tags and bool(tags & {"path", "polyline", "polygon", "line", "rect"}),
            "missing SVG text or vector graphics")
    require("image" not in tags, "unexpected raster image in SVG")


def check_pdf(path, width, height):
    pdf = path.read_bytes()
    require(pdf.startswith(b"%PDF-") and pdf.rstrip().endswith(b"%%EOF"),
            "invalid PDF header or trailer")
    info = subprocess.run(
        ["pdfinfo", str(path.resolve())], capture_output=True, text=True, timeout=10,
        env={**os.environ, "LC_ALL": "C"},
    )
    require(info.returncode == 0, info.stderr.strip() or "could not read PDF")
    require(bool(re.search(r"^Pages:[ \t]+1$", info.stdout, re.MULTILINE)),
            "expected a single PDF page")
    page_size = re.search(r"^Page size:[ \t]+([0-9.]+) x ([0-9.]+) pts", info.stdout, re.MULTILINE)
    require(page_size is not None, "missing PDF page size")
    actual = tuple(float(value) for value in page_size.groups())
    expected = (width * POINTS_PER_MM, height * POINTS_PER_MM)
    require(all(abs(value - target) < 1.01 for value, target in zip(actual, expected)),
            f"expected PDF size {expected} pt, found {actual}")



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(__file__).resolve().parents[1] / "build" / "figures"
    parser.add_argument("directory", nargs="?", type=Path, default=default)
    args = parser.parse_args()
    if shutil.which("pdfinfo") is None:
        parser.exit(1, "Install Poppler (pdfinfo) to check PDF figures.\n")
    checks = {"pdf": check_pdf, "svg": check_svg, "png": check_png}
    sizes = expected_sizes()
    failures = []
    for stem, (width, height) in sizes.items():
        for extension, check in checks.items():
            filename = f"{stem}.{extension}"
            try:
                check(args.directory / filename, width, height)
            except (OSError, ValueError, KeyError, struct.error, ET.ParseError,
                    subprocess.TimeoutExpired) as error:
                failures.append(f"{filename}: {error}")
    if failures:
        parser.exit(1, "\n".join(failures) + "\n")
    print(f"Validated {len(PLATES)} figure plates and {len(sizes) - len(PLATES)} panels "
          f"in PDF, SVG and {DPI}-dpi PNG ({len(sizes) * len(checks)} files).")


if __name__ == "__main__":
    main()
