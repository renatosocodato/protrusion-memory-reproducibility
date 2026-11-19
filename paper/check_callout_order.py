#!/usr/bin/env python3
"""Check first-use order and figure-display links for this Markdown manuscript."""
from pathlib import Path
import re
import sys

p = Path(__file__).with_name('MANUSCRIPT_WORKING_DRAFT_EN.md')
s = p.read_text(encoding='utf8')
pattern = re.compile(r'\b(Supplementary Figure|Supplementary Table|Supplementary Section|Figure|Table|Equation)\s+(\d+)([a-z]?)\b(?:, panel ([a-z])\b)?')
seen = {}
errors = []
for match in pattern.finditer(s):
    series, num, panel = match.group(1), int(match.group(2)), match.group(3)
    series_seen = seen.setdefault(series, [])
    panels = (panel, match.group(4)) if match.group(4) else (panel,)
    for letter in panels:
        item = (num, letter)
        if item not in series_seen:
            if series_seen and item < series_seen[-1]:
                errors.append(f'{series} {num}{letter} first appears after {series_seen[-1]}')
            series_seen.append(item)
# Grouped numbers and ranges require explicit expanded-order checking.
if re.search(r'\b(?:Figure|Table|Equation)s?\s+\d+[a-z]?\s*(?:[–-]|,|and)\s*\d',s):
    errors.append('Grouped/range call requires explicit expanded-order support')
figure_captions = [int(x) for x in re.findall(r'\*\*Figure\s+(\d+)\.',s)]
if figure_captions != sorted(set(figure_captions)):
    errors.append(f'Figure captions have wrong order: {figure_captions}')
if figure_captions != sorted({n for n, _ in seen.get('Figure', [])}):
    errors.append('Figure captions and substantive mentions differ')
for link in re.findall(r'!\[[^]]*\]\(([^)]+)\)',s):
    if not (p.parent / link).resolve().is_file(): errors.append('Missing image: '+link)
if errors:
    print('\n'.join(errors));sys.exit(1)
print('First-use series:',seen,'; figure display order:',figure_captions,'; image links resolve.')
