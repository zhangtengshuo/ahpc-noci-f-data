# INDO-dimer performance

Figs. 4/5; SI S8/S9/S11.

## Software-release run examples

[release_version_run_examples/](release_version_run_examples/README.md) contains the inputs and outputs for calculations with NOCI.jl 26.10.0, including monomer states, CommonOrb bases and NOCI matrices and energies.
The existing `server/`, `laptop/`, `inputs/`, `monomer/` and `environment/` contain the original datasets from earlier development versions; the remaining sections of this README describe those earlier data.

## Original development datasets

The MEBF order is S0S0, S1S0, S0S1, S1S1, T1T1, D-D+, D+D-.
In run names, `0g8t` denotes zero GPUs and eight physical CPU cores; `8g8t` denotes eight GPUs and eight physical CPU cores.
`server/` contains numerical outputs and timing logs; `laptop/` contains the laptop input and outputs.
The logs contain full-process timing and peak memory; `telemetry.phases` in JSON files contains stage timings.
`inputs/` contains the NOCI inputs, fragment wave functions, geometry, CommonOrb, one-electron integrals and SI Molden file.
`monomer/` contains the five state-specific RasOrb files.
The two-electron factor bundle is not included.

Recalculate reuse counts with `julia --project=/path/to/NOCI.jl ../reuse_counts.jl inputs/INDO_0g1t.yaml OUTPUT.csv`.
The 28 INDO rows are in [occupation_string_pairs.csv](../04_pentacene/occupation_string_pairs.csv), selected by `system=INDO`.
The published GronOR reference is Table 6 of Straatsma et al., JCTC 2022, DOI `10.1021/acs.jctc.2c00266`.
`environment/` supplies the Julia dependency files; see [SOFTWARE.md](../SOFTWARE.md).

## File list

| File | Contents / use |
|---|---|
| [environment/Manifest.toml](environment/Manifest.toml) | Locked Julia dependency versions. |
| [environment/Project.toml](environment/Project.toml) | Julia package specification. |
| [inputs/B1_INDOAB_GS0S0.xyz](inputs/B1_INDOAB_GS0S0.xyz) | Dimer geometry, in angstroms. |
| [inputs/B2_COMMONORB_INDOAB_GS0S0_CB1e-5](inputs/B2_COMMONORB_INDOAB_GS0S0_CB1e-5) | Aggregate CommonOrb basis. |
| [inputs/B2_INDOA_D+@GS0.det](inputs/B2_INDOA_D+@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOA_D+@GS0.vec](inputs/B2_INDOA_D+@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOA_D-@GS0.det](inputs/B2_INDOA_D-@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOA_D-@GS0.vec](inputs/B2_INDOA_D-@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOA_S0@GS0.det](inputs/B2_INDOA_S0@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOA_S0@GS0.vec](inputs/B2_INDOA_S0@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOA_S1@GS0.det](inputs/B2_INDOA_S1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOA_S1@GS0.vec](inputs/B2_INDOA_S1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOA_T1@GS0.det](inputs/B2_INDOA_T1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOA_T1@GS0.vec](inputs/B2_INDOA_T1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOB_D+@GS0.det](inputs/B2_INDOB_D+@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOB_D+@GS0.vec](inputs/B2_INDOB_D+@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOB_D-@GS0.det](inputs/B2_INDOB_D-@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOB_D-@GS0.vec](inputs/B2_INDOB_D-@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOB_S0@GS0.det](inputs/B2_INDOB_S0@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOB_S0@GS0.vec](inputs/B2_INDOB_S0@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOB_S1@GS0.det](inputs/B2_INDOB_S1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOB_S1@GS0.vec](inputs/B2_INDOB_S1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_INDOB_T1@GS0.det](inputs/B2_INDOB_T1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_INDOB_T1@GS0.vec](inputs/B2_INDOB_T1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/INDOAB.one](inputs/INDOAB.one) | One-electron integrals and nuclear repulsion. |
| [inputs/INDO_0g1t.yaml](inputs/INDO_0g1t.yaml) | NOCI input; CPU-only settings are identical for all six CPU allocations. |
| [inputs/INDO_1g1t.yaml](inputs/INDO_1g1t.yaml) | NOCI input for this CPU/GPU allocation. |
| [inputs/INDO_2g2t.yaml](inputs/INDO_2g2t.yaml) | NOCI input for this CPU/GPU allocation. |
| [inputs/INDO_4g4t.yaml](inputs/INDO_4g4t.yaml) | NOCI input for this CPU/GPU allocation. |
| [inputs/INDO_8g8t.yaml](inputs/INDO_8g8t.yaml) | NOCI input for this CPU/GPU allocation. |
| [inputs/INDO_monomer_S0_CAS88_vdZP.rasscf.molden](inputs/INDO_monomer_S0_CAS88_vdZP.rasscf.molden) | Monomer orbital source for the SI active-space figure. |
| [laptop/08_INDOAB_GS0S0_CB1e-5_AHPC_CHO_12CPU.day](laptop/08_INDOAB_GS0S0_CB1e-5_AHPC_CHO_12CPU.day) | Laptop event log, including solve and completion times. |
| [laptop/08_INDOAB_GS0S0_CB1e-5_AHPC_CHO_12CPU.json](laptop/08_INDOAB_GS0S0_CB1e-5_AHPC_CHO_12CPU.json) | Laptop numerical output. |
| [laptop/input.yaml](laptop/input.yaml) | Laptop input. |
| [monomer/03_MolX_D+@GS0_C1_desym.RasOrb](monomer/03_MolX_D+@GS0_C1_desym.RasOrb) | D+ monomer RasOrb at GS0. |
| [monomer/03_MolX_D-@GS0_C1_desym.RasOrb](monomer/03_MolX_D-@GS0_C1_desym.RasOrb) | D- monomer RasOrb at GS0. |
| [monomer/03_MolX_S0@GS0_C1_desym.RasOrb](monomer/03_MolX_S0@GS0_C1_desym.RasOrb) | S0 monomer RasOrb at GS0. |
| [monomer/03_MolX_S1@GS0_C1_desym.RasOrb](monomer/03_MolX_S1@GS0_C1_desym.RasOrb) | S1 monomer RasOrb at GS0. |
| [monomer/03_MolX_T1@GS0_C1_desym.RasOrb](monomer/03_MolX_T1@GS0_C1_desym.RasOrb) | T1 monomer RasOrb at GS0. |
| [server/INDO_0g16t.json](server/INDO_0g16t.json) | Matrices, energies, determinant counts and stage timings for INDO_0g16t. |
| [server/INDO_0g16t_1327.out](server/INDO_0g16t_1327.out) | Event log and GNU-time full-process timing/memory for INDO_0g16t. |
| [server/INDO_0g1t.json](server/INDO_0g1t.json) | Matrices, energies, determinant counts and stage timings for INDO_0g1t. |
| [server/INDO_0g1t_1331.out](server/INDO_0g1t_1331.out) | Event log and GNU-time full-process timing/memory for INDO_0g1t. |
| [server/INDO_0g2t.json](server/INDO_0g2t.json) | Matrices, energies, determinant counts and stage timings for INDO_0g2t. |
| [server/INDO_0g2t_1330.out](server/INDO_0g2t_1330.out) | Event log and GNU-time full-process timing/memory for INDO_0g2t. |
| [server/INDO_0g32t.json](server/INDO_0g32t.json) | Matrices, energies, determinant counts and stage timings for INDO_0g32t. |
| [server/INDO_0g32t_1326.out](server/INDO_0g32t_1326.out) | Event log and GNU-time full-process timing/memory for INDO_0g32t. |
| [server/INDO_0g4t.json](server/INDO_0g4t.json) | Matrices, energies, determinant counts and stage timings for INDO_0g4t. |
| [server/INDO_0g4t_1329.out](server/INDO_0g4t_1329.out) | Event log and GNU-time full-process timing/memory for INDO_0g4t. |
| [server/INDO_0g8t.json](server/INDO_0g8t.json) | Matrices, energies, determinant counts and stage timings for INDO_0g8t. |
| [server/INDO_0g8t_1328.out](server/INDO_0g8t_1328.out) | Event log and GNU-time full-process timing/memory for INDO_0g8t. |
| [server/INDO_1g1t.json](server/INDO_1g1t.json) | Matrices, energies, determinant counts and stage timings for INDO_1g1t. |
| [server/INDO_1g1t_1325.out](server/INDO_1g1t_1325.out) | Event log and GNU-time full-process timing/memory for INDO_1g1t. |
| [server/INDO_2g2t.json](server/INDO_2g2t.json) | Matrices, energies, determinant counts and stage timings for INDO_2g2t. |
| [server/INDO_2g2t_1324.out](server/INDO_2g2t_1324.out) | Event log and GNU-time full-process timing/memory for INDO_2g2t. |
| [server/INDO_4g4t.json](server/INDO_4g4t.json) | Matrices, energies, determinant counts and stage timings for INDO_4g4t. |
| [server/INDO_4g4t_1323.out](server/INDO_4g4t_1323.out) | Event log and GNU-time full-process timing/memory for INDO_4g4t. |
| [server/INDO_8g8t.json](server/INDO_8g8t.json) | Matrices, energies, determinant counts and stage timings for INDO_8g8t. |
| [server/INDO_8g8t_1322.out](server/INDO_8g8t_1322.out) | Event log and GNU-time full-process timing/memory for INDO_8g8t. |
