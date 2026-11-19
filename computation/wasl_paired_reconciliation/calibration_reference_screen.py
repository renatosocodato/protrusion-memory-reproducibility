#!/usr/bin/env python3
"""Exploratory descriptive reference screen; no inferred cell replication."""
import hashlib
import json
from pathlib import Path
from reproduction_io import amiri_workbook, output_directory

import numpy as np
from openpyxl import load_workbook

HERE = output_directory(Path(__file__))
WORKBOOK = amiri_workbook()


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    w = load_workbook(WORKBOOK, read_only=True, data_only=True)
    sheet = w['SFig 8']
    headers = [sheet.cell(2, i).value for i in range(2, 8)]
    expected = ['CTS', 'TS', 'TS+calyc', 'TS + bleb', 'Lysed TS', 'Lysed CTS']
    if headers != expected:
        raise ValueError(f'Source headers differ from frozen protocol: {headers}')
    conditions = {}
    for column, label in enumerate(headers, start=2):
        data = []
        for row in sheet.iter_rows(min_row=3, min_col=column, max_col=column):
            value = row[0].value
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                data.append(float(value))
        x = np.asarray(data)
        if not len(x) or not np.isfinite(x).all() or np.min(x)<0 or np.max(x)>1:
            raise ValueError(f'Unexpected FRET efficiency range for {label}')
        conditions[label] = {'n_values': len(x), 'min': float(np.min(x)), 'median': float(np.median(x)),
                             'q25': float(np.quantile(x,.25)), 'q75': float(np.quantile(x,.75)),
                             'max': float(np.max(x))}
    w.close()
    control = conditions['CTS']['median']
    for label in ('Lysed TS','Lysed CTS'):
        conditions[label]['median_difference_vs_CTS_percentage_points'] = 100*(conditions[label]['median']-control)
    out = {'source': 'Amiri et al. 2023 SFig 8 workbook; U2OS mechanical sensor; no N-WASP',
           'workbook_sha256': sha(WORKBOOK),
           'row_unit': 'Unknown; measurements may include repeated time points per cell',
           'conditions': conditions,
           'interpretation_limit': 'Distributions summarize reported values without cell/donor independence or FRET-to-force calibration.'}
    (HERE/'calibration_reference_screen.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
