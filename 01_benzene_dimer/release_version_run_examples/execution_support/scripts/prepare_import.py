#!/usr/bin/env python3
"""Render official NOCI run inputs from 00 definitions and imported artifacts."""
import copy
import json
from pathlib import Path
import sys
import yaml

root = Path(sys.argv[1]).resolve()
runs = json.loads((root / "00_plan/runs.json").read_text())
for run in runs:
    case = root / run["case"]
    raw = yaml.safe_load((case / f"00_system_define/{run['variant']}.yaml").read_text())
    reuse = yaml.safe_load((case / "00_system_define/reuse.yaml").read_text())
    table = raw["mebf_table"]
    mebfs = []
    used = {}
    for column in range(1, len(table[0])):
        labels = []
        for row in table:
            item = reuse["bindings"][f"{row[0]}/{row[column]}"]
            labels.append(item["label"])
            used[item["label"]] = item
        mebf = dict(label=reuse["mebf_labels"][column-1], states=labels)
        coupling = reuse["coupling_multiplicities"][column-1]
        if coupling is not None:
            mebf["coupling_multiplicities"] = coupling
        mebfs.append(mebf)
    gpu = raw["runtime"].get("use_gpu", False)
    runtime = dict(use_gpu=gpu, scratch_bytes=raw["runtime"]["scratch_mb"]*1024**2,
                   max_host_bytes=(48 if gpu else 24 if raw["runtime"]["threads"] == 8 else 12)*1024**3,
                   collect_allocations=False)
    if gpu:
        runtime.update(gpu_devices=raw["runtime"]["gpu_devices"], max_device_bytes=10*1024**3)
    if run["name"].startswith(("benzene", "indo")):
        runtime["max_active_occupied"] = 12
    noci = dict(algorithm=dict(name="ahpc-noci", branch="spin_string", support_lookup="adaptive"),
                system=dict(label=raw["project"], spin_multiplicity=raw["system"]["spin_multiplicity"],
                            fragment_states=list(used.values()), mebfs=mebfs),
                integrals=copy.deepcopy(reuse["integrals_packed" if run["storage"] == "packed" else "integrals_factor"]),
                screening=dict(raw["screening"], tau_singularity=reuse["tau_singularity"]), runtime=runtime)
    output = case / "05_noci" / run["variant"]
    output.mkdir(parents=True, exist_ok=True)
    (output / "results").mkdir(exist_ok=True)
    (output / "input.yaml").write_text(yaml.safe_dump(noci, sort_keys=False))
    receipt = dict(definition=str(case / f"00_system_define/{run['variant']}.yaml"),
                   reuse_map=str(case / "00_system_define/reuse.yaml"),
                   julia_threads=raw["runtime"]["threads"], blas_threads=raw["runtime"]["blas_threads"],
                   commonorb_threshold=raw["common_orbitals"]["threshold"],
                   imported_integral_threshold=reuse["imported_integral_threshold"],
                   preparation="historical products imported; no fabricated stage completion records")
    (output / "preparation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(output / "input.yaml")
