# Data for AHPC-NOCI-F

Data supporting “Efficient Nonorthogonal Configuration Interaction via Active-Space Hamiltonian Precontraction and Polynomial Determinant-Derivative Contraction”.

| Directory | Contents | Paper location |
|---|---|---|
| [01_benzene_dimer](01_benzene_dimer/README.md) | Benzene-dimer validation, CommonOrb and monomer RasOrb | Main-text Section 4.1; SI S6/S11 |
| [02_pddc](02_pddc/README.md) | PDDC accuracy, timing and storage | Fig. 3; SI S7 |
| [03_indo_dimer](03_indo_dimer/README.md) | INDO-dimer CPU/GPU performance, CommonOrb and monomer RasOrb | Figs. 4/5; SI S8/S9/S11 |
| [04_pentacene](04_pentacene/README.md) | Energies, matrices, couplings, reuse counts and orbitals | Figs. 5/6; SI S9/S10/S11 |
| [05_ethene_example](05_ethene_example/README.md) | AHPC numerical example and complete occupation-pattern reuse enumeration | SI S2/S3; Fig. S1 |

## Release-version runs and original development data

The following directories contain calculation inputs and outputs obtained with NOCI.jl 26.10.0:

- [Benzene release-version run examples](01_benzene_dimer/release_version_run_examples/README.md): two single-core AHPC runs, using packed integrals and factors from the same packed data.
- [INDO release-version run examples](03_indo_dimer/release_version_run_examples/README.md): 8-core CPU and 8-core/8-GPU runs.
- [Pentacene release-version run examples](04_pentacene/release_version_run_examples/README.md): dimer and trimer, each using 8 physical CPU cores and 8 GPUs.

The other existing calculation folders retain the original datasets from earlier development versions.
Stages 01–04 contain monomer states, geometries, CommonOrb bases and integral records; stage 05 contains the release-version NOCI inputs, matrices and energies.
CommonOrb bases are included for all molecular examples.

Each directory README lists its files and calculation commands.
[SOFTWARE.md](SOFTWARE.md) describes the software and shared scripts.
