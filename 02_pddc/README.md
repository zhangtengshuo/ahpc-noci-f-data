# PDDC accuracy, timing and storage

Fig. 3; SI S7.

`accuracy/accuracy.csv` contains eight controlled six-orbital cases and twenty diagonal analytic cases.
Here $q$ is the number of occupied active-space orbitals in one spin channel.
The columns `c0`, `c1` and `r2` refer to the determinant, first-cofactor matrix and contracted second-order scalar.
`r2_generic` denotes the ordered-index contraction; `r2_hessian` denotes the compressed pair-index evaluation.
These algebraic test quantities are dimensionless.
`timing/timing.csv` contains five measurements at each of seven dimensions; `seconds_per_task` is in seconds.
`timing/memory.csv` reports storage in bytes.

Run accuracy with `julia --project=/path/to/NOCI.jl -t 1 accuracy/section42_pddc_experiment.jl OUTPUT_DIRECTORY accuracy`.
Run timing/storage with `julia --project=/path/to/NOCI.jl -t 1 timing/section42_pddc_experiment.jl OUTPUT_DIRECTORY timing`.
Set `OPENBLAS_NUM_THREADS=1` for both commands.
`environment/` supplies the Julia dependency files; see [SOFTWARE.md](../SOFTWARE.md).

## File list

| File | Contents / use |
|---|---|
| [accuracy/accuracy.csv](accuracy/accuracy.csv) | 28-case accuracy data for Fig. 3a and SI S7. |
| [accuracy/run.json](accuracy/run.json) | Experiment/software/hardware metadata. |
| [accuracy/section42_pddc_experiment.jl](accuracy/section42_pddc_experiment.jl) | Experiment script: matrix construction, reference calculations and measurement procedure. |
| [environment/Manifest.toml](environment/Manifest.toml) | Locked Julia dependency versions. |
| [environment/Project.toml](environment/Project.toml) | Julia package specification. |
| [timing/memory.csv](timing/memory.csv) | Measured polynomial workspace and coefficient-tensor storage. |
| [timing/run.json](timing/run.json) | Experiment/software/hardware metadata. |
| [timing/section42_pddc_experiment.jl](timing/section42_pddc_experiment.jl) | Experiment script: matrix construction, reference calculations and measurement procedure. |
| [timing/timing.csv](timing/timing.csv) | All 35 single-core timing measurements for Fig. 3b and SI S7. |
