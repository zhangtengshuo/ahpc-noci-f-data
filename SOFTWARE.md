# Software and running the calculations

The original development-version molecular and PDDC scripts use Julia 1.12.6 and [NOCI.jl](https://github.com/zhangtengshuo/NOCI.jl).
Use the corresponding source revision and the supplied dependency files with that checkout.

| Calculations | NOCI.jl source revision | Dependency files |
|---|---|---|
| Benzene, PDDC and ethylene scripts | `fe4ae809257522f6f87b640911c7e0dcb9f0c5fd` | [02_pddc/environment](02_pddc/environment/Project.toml) |
| INDO | `3bfe4fb9d0d7488055564356daf29a862bcf208d` | [03_indo_dimer/environment](03_indo_dimer/environment/Project.toml) |
| Pentacene | `46f61d4` | [04_pentacene/environment](04_pentacene/environment/Project.toml) |

Install the Julia dependencies once with `julia --project=/path/to/NOCI.jl -e 'using Pkg; Pkg.instantiate()'`.
Use `--project=/path/to/NOCI.jl` when running the Julia scripts.
Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`; select Julia threads with `-t N`.
For CPU timing, bind those threads to distinct physical CPU cores.
GPU runs additionally require the devices selected by the NOCI input.

| Shared script | Usage |
|---|---|
| [run_noci.jl](run_noci.jl) | `julia --project=/path/to/NOCI.jl -t N run_noci.jl INPUT.yaml OUTPUT_DIRECTORY LABEL` |
| [reuse_counts.jl](reuse_counts.jl) | `julia --project=/path/to/NOCI.jl -t N reuse_counts.jl INPUT.yaml OUTPUT.csv` |
| [analyze_matrices.py](analyze_matrices.py) | `python analyze_matrices.py RESULT.json OUTPUT_DIRECTORY` or use a benzene result `.txt` |
| [dataset_input.jl](dataset_input.jl) | Input resolver and coupling expressions used by the Julia scripts |

The input-path resolver locates supplied files by basename within each case and leaves the archived YAML files unchanged.
The count script uses fragment wave functions and screening settings without loading two-electron integrals.
The ethylene case includes the integrals needed to run its NOCI calculation.
The benzene, INDO and pentacene cases retain their numerical outputs and count inputs; their two-electron datasets are not distributed.

`analyze_matrices.py` requires Python 3 and NumPy and exports energy and pairwise-coupling CSV files.
It uses the normalized matrices and the two-state symmetric-orthogonalization expression from the manuscript.
For pentacene it processes both the unshifted and CASPT2-shifted Hamiltonians.
State phases can be set with repeated `--flip-state LABEL` options.
`05_ethene_example/reuse/enumerate_reuse.py` requires only Python 3.

## Software-release run examples

The release-version examples use NOCI.jl 26.10.0, source revision `d3c44af10db6d3ef483544f803c5c6a24774eb4e`, Julia 1.12.6 and OpenMolcas4NOCI 26.06.2.
The calculations use the `noci run` entry point and the same monomer wave functions and integral data as the earlier development calculations.
The input files specify the MEBF states, screening thresholds and CPU/GPU settings; the outputs contain NOCI matrices, energies and CASPT2 corrections.
The dependency tables and shared-script instructions above describe the original development datasets.
