#!/bin/bash
set -euo pipefail
case_dir="/scratch/alpha/AHPC/$1"
variant="$2"
threads="$3"
gpu="$4"
convert="$5"
source /home/alpha/bin/NOCI.jl/v26.10.0/bin/noci-environment.sh
export JULIA_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MOLCAS_NPROCS=1
export JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_LOAD_PATH='@:@stdlib'
project="$NOCI_ROOT"
if [[ "$gpu" == true ]]; then project="$NOCI_ROOT/optional/gpu"; fi
export JULIA_PROJECT="$project" NOCI_JULIA_PROJECT="$project"
output="$case_dir/05_noci/$variant"
scratch_root=""
cleanup() {
    result=$?
    trap - EXIT
    printf '%s\n' "$result" > "$output/exit_status.txt"
    if [[ -n "$scratch_root" ]]; then rm -rf -- "$scratch_root"; fi
    exit "$result"
}
trap cleanup EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
cd "$output"
python3 - "$output" <<'PY'
import datetime,json,os,pathlib,sys
p=pathlib.Path(sys.argv[1])
receipt=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    slurm_job_id=os.environ.get('SLURM_JOB_ID'), allocated_cpus=sorted(os.sched_getaffinity(0)),
    julia_threads=os.environ['JULIA_NUM_THREADS'], blas_threads=os.environ['OPENBLAS_NUM_THREADS'],
    cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'), project=os.environ['JULIA_PROJECT'],
    installation='/home/alpha/bin/NOCI.jl/v26.10.0',
    installed_source_commit='d3c44af10db6d3ef483544f803c5c6a24774eb4e')
(p/'execution.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2),flush=True)
PY
if [[ "$convert" == true ]]; then
    scratch_root=$(mktemp -d "/dev/shm/AHPC_${SLURM_JOB_ID}_XXXXXX")
    python3 - "$case_dir" "$scratch_root" "$SLURM_JOB_ID" <<'PY'
import pathlib,sys,yaml
p=pathlib.Path(sys.argv[1])/'04_integrals'
x=yaml.safe_load((p/'seward_to_factors.yaml').read_text())
x['seward']['scratch_parent']=sys.argv[2]
(p/('seward_to_factors_job_'+sys.argv[3]+'.yaml')).write_text(yaml.safe_dump(x,sort_keys=False))
PY
    # The installed GPU environment omits a direct YAML dependency; use the installed CPU project as its explicit secondary environment for this public conversion tool.
    JULIA_LOAD_PATH="@:$NOCI_ROOT:@stdlib" /usr/bin/time -v -o "$case_dir/04_integrals/conversion.time" \
        "$JULIA_BIN" --startup-file=no --project="$project" --threads="$threads" --heap-size-hint=16G \
        "$NOCI_ROOT/tools/seward_to_factor_bundle.jl" \
        "$case_dir/04_integrals/seward_to_factors_job_${SLURM_JOB_ID}.yaml" \
        "$case_dir/04_integrals/TriPAB_cho_factor_bundle.jld2" \
        > "$case_dir/04_integrals/conversion.log" 2>&1
fi
/usr/bin/time -v -o "$output/run.time" \
    "$JULIA_BIN" --startup-file=no --project="$project" --threads="$threads" --heap-size-hint=16G \
    "$NOCI_ROOT/bin/noci" run "$output/input.yaml" "$output/results" result \
    "$case_dir/00_system_define/geometry/aggregate.xyz" \
    > "$output/run.log" 2>&1
date -Is > "$output/completed_at.txt"
