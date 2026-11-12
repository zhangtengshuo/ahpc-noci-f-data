#!/usr/bin/env python3
"""Collect operational receipts and output paths without reading scientific results."""
import datetime
import json
from pathlib import Path
import subprocess
import sys
import tarfile

root = Path(sys.argv[1]).resolve()
runs = json.loads((root/'00_plan/runs.json').read_text())
submissions = json.loads((root/'00_plan/submissions.json').read_text())
jobs = {x['name']: x['job_id'] for x in submissions}
records = []
for run in runs:
    work = root/run['case']/'05_noci'/run['variant']
    status = work/'exit_status.txt'
    execution = work/'execution.json'
    timing = work/'run.time'
    metrics = {}
    if timing.exists():
        for line in timing.read_text().splitlines():
            for key in ('Elapsed (wall clock) time','Maximum resident set size','Exit status'):
                if line.strip().startswith(key):
                    metrics[key] = line.strip()
    records.append(dict(name=run['name'],job_id=jobs[run['name']],directory=str(work),
                        exit_status=int(status.read_text()) if status.exists() else None,
                        execution=json.loads(execution.read_text()) if execution.exists() else None,
                        timing=metrics,
                        output_files=[dict(name=p.name,bytes=p.stat().st_size) for p in sorted((work/'results').glob('*')) if p.is_file()]))
accounting = subprocess.run(['sacct','-j',','.join(jobs.values()),'--format=JobID,JobName%32,State,Elapsed,ExitCode','-X','-n'],capture_output=True,text=True,check=True).stdout
receipt = dict(collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               operational_records=records,scheduler_accounting=accounting,
               scientific_comparison_performed=False)
(root/'00_plan/completion.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(accounting)
with tarfile.open(root/'00_plan/inputs_outputs_logs.tar.gz','w:gz') as tar:
    files = [p for p in (root/'00_plan').glob('*') if p.is_file() and p.suffix in ('.json','.md')]
    files.extend((root/'logs').glob('*.log'))
    for run in runs:
        case = root/run['case']
        files.extend(p for p in (case/'00_system_define').glob('*.yaml') if p.is_file())
        work = case/'05_noci'/run['variant']
        files.extend(p for p in work.iterdir() if p.is_file())
        files.extend(p for p in (work/'results').iterdir() if p.is_file() and p.suffix in ('.out','.json','.sec'))
        if run['conversion']:
            files.extend(p for p in (case/'04_integrals').iterdir() if p.is_file() and not p.is_symlink() and p.suffix in ('.log','.time','.yaml'))
    for p in dict.fromkeys(files):
        tar.add(p,arcname=str(p.relative_to(root)),recursive=False)
print('Operational receipt and text input/output archive written; no scientific result comparison.')
