#!/usr/bin/env python3
"""Materialize explicit reuse links and retain printed settings as provenance."""
import json
from pathlib import Path
import re
import sys

root = Path(sys.argv[1]).resolve()
links = json.loads((root / '00_plan/links.json').read_text())
for item in links:
    target = root / item['destination']
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.is_symlink():
        target.symlink_to(item['source'])
(root/'logs').mkdir(exist_ok=True)
evidence = []
for item in links:
    if not item['destination'].endswith(('.out','.input')):
        continue
    source = Path(item['source'])
    if not source.is_file():
        continue
    lines = source.read_text(errors='replace').splitlines()
    selected = set()
    for i,line in enumerate(lines):
        if re.search(r'CD Threshold|acCD|Cholesky decomposition threshold|CommonOrb|Threshold|aCD auxiliary',line,re.I):
            selected.update(range(max(0,i-1),min(len(lines),i+2)))
    if selected:
        evidence.append(dict(source=str(source), excerpts=[dict(line=i+1,text=lines[i]) for i in sorted(selected)]))
(root/'00_plan/printed_settings.json').write_text(json.dumps(evidence,indent=2)+'\n')
print('Historical reuse links and printed-setting excerpts prepared; no calculation executed.')
