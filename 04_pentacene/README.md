# Pentacene dimer and trimer

Figs. 5/6; SI S9/S10/S11.

## Software-release run examples

[release_version_run_examples/](release_version_run_examples/README.md) contains the inputs and outputs from the installed NOCI.jl 26.10.0 release calculations, organized as complete 00–05 workflows.
The existing `dimer/`, `trimer/`, `monomer/` and `environment/` directories and the original reuse-count table contain the original datasets from earlier development versions; the remaining sections of this README describe those earlier data.

## Original development datasets

`dimer/` and `trimer/` contain geometries, fragment wave functions, NOCI inputs and numerical outputs.
The dimer and trimer CommonOrb bases are included.
`monomer/` contains the state-specific RasOrb files, the two initial RasOrb files, the CASSCF/CASPT2 inputs/outputs, monomer geometries and SI Molden file.
Two-electron integral datasets are not included.

`*_normalized_hamiltonian_Eh.csv` and `*_normalized_overlap.csv` contain the unshifted SI matrices.
`*_caspt2_shifted_hamiltonian_Eh.csv` contains the shifted Hamiltonian used for Fig. 6.
Matrix labels identify the MEBFs; Hamiltonian elements and energies are in Hartree, and overlaps are dimensionless.
The result JSON contains raw matrices, NOCI eigenvalues and CASPT2 corrections in `baseline_caspt2.mebf_shifts`.
The trimer data retain the fragment-state phases used in the SI; Fig. 6 flips both CT-state phases.

`occupation_string_pairs.csv` contains the Fig. 5/SI S9 reuse data: 28 INDO, 21 dimer and 21 trimer rows.
Recalculate the dimer counts with `julia --project=/path/to/NOCI.jl ../reuse_counts.jl dimer/dimer.yaml OUTPUT.csv`; use `trimer/trimer.yaml` for the trimer.
Analyze a result with `python ../analyze_matrices.py RESULT.json OUTPUT_DIRECTORY`.
For the Fig. 6 trimer phase convention, add `--flip-state D+D-S0 --flip-state D-D+S0`.
`environment/` supplies the Julia dependency files; see [SOFTWARE.md](../SOFTWARE.md).

## File list

| File | Contents / use |
|---|---|
| [dimer/B1_TriPAB_GS1S0.xyz](dimer/B1_TriPAB_GS1S0.xyz) | Aggregate geometry, in angstroms. |
| [dimer/B2_COMMONORB_TriPAB_GS1S0_CB1e-5](dimer/B2_COMMONORB_TriPAB_GS1S0_CB1e-5) | Aggregate CommonOrb basis. |
| [dimer/B2_TriPA_D+@GS1.det](dimer/B2_TriPA_D+@GS1.det) | Fragment determinant expansion. |
| [dimer/B2_TriPA_D+@GS1.vec](dimer/B2_TriPA_D+@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPA_D-@GS1.det](dimer/B2_TriPA_D-@GS1.det) | Fragment determinant expansion. |
| [dimer/B2_TriPA_D-@GS1.vec](dimer/B2_TriPA_D-@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPA_S0@GS1.det](dimer/B2_TriPA_S0@GS1.det) | Fragment determinant expansion. |
| [dimer/B2_TriPA_S0@GS1.vec](dimer/B2_TriPA_S0@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPA_S1@GS1.det](dimer/B2_TriPA_S1@GS1.det) | Fragment determinant expansion. |
| [dimer/B2_TriPA_S1@GS1.vec](dimer/B2_TriPA_S1@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPA_T1@GS1.det](dimer/B2_TriPA_T1@GS1.det) | Fragment determinant expansion. |
| [dimer/B2_TriPA_T1@GS1.vec](dimer/B2_TriPA_T1@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPB_D+@GS0.det](dimer/B2_TriPB_D+@GS0.det) | Fragment determinant expansion. |
| [dimer/B2_TriPB_D+@GS0.vec](dimer/B2_TriPB_D+@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPB_D-@GS0.det](dimer/B2_TriPB_D-@GS0.det) | Fragment determinant expansion. |
| [dimer/B2_TriPB_D-@GS0.vec](dimer/B2_TriPB_D-@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPB_S0@GS0.det](dimer/B2_TriPB_S0@GS0.det) | Fragment determinant expansion. |
| [dimer/B2_TriPB_S0@GS0.vec](dimer/B2_TriPB_S0@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPB_S1@GS0.det](dimer/B2_TriPB_S1@GS0.det) | Fragment determinant expansion. |
| [dimer/B2_TriPB_S1@GS0.vec](dimer/B2_TriPB_S1@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/B2_TriPB_T1@GS0.det](dimer/B2_TriPB_T1@GS0.det) | Fragment determinant expansion. |
| [dimer/B2_TriPB_T1@GS0.vec](dimer/B2_TriPB_T1@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [dimer/TriPAB.one](dimer/TriPAB.one) | One-electron integrals and nuclear repulsion. |
| [dimer/TriPAB_GS1S0_CB1e-5_CHO_Low_tau1e-5.json](dimer/TriPAB_GS1S0_CB1e-5_CHO_Low_tau1e-5.json) | Raw/normalized matrices, energies, CASPT2 corrections and timings. |
| [dimer/TriPAB_GS1S0_CB1e-5_CHO_Low_tau1e-5.occupation_statistics.json](dimer/TriPAB_GS1S0_CB1e-5_CHO_Low_tau1e-5.occupation_statistics.json) | 21-pair occupation-string statistics for Fig. 5 and SI S9. |
| [dimer/TriPAB_GS1S0_CB1e-5_CHO_Low_tau1e-5.sec](dimer/TriPAB_GS1S0_CB1e-5_CHO_Low_tau1e-5.sec) | Calculation event log. |
| [dimer/dimer.yaml](dimer/dimer.yaml) | NOCI input, with both screening thresholds equal to 1e-5. |
| [dimer/dimer_caspt2_shifted_hamiltonian_Eh.csv](dimer/dimer_caspt2_shifted_hamiltonian_Eh.csv) | CASPT2-shifted normalized Hamiltonian in Hartree; Fig. 6. |
| [dimer/dimer_normalized_hamiltonian_Eh.csv](dimer/dimer_normalized_hamiltonian_Eh.csv) | Unshifted normalized Hamiltonian in Hartree; SI S10. |
| [dimer/dimer_normalized_overlap.csv](dimer/dimer_normalized_overlap.csv) | Normalized overlap matrix; SI S10. |
| [environment/Manifest.toml](environment/Manifest.toml) | Locked Julia dependency versions. |
| [environment/Project.toml](environment/Project.toml) | Julia package specification. |
| [monomer/01_TriP_CAS66_S0_seed.RasOrb](monomer/01_TriP_CAS66_S0_seed.RasOrb) | S0 initial RasOrb referenced by the monomer inputs. |
| [monomer/01_TriP_CAS66_S1_seed.RasOrb](monomer/01_TriP_CAS66_S1_seed.RasOrb) | S1 initial RasOrb referenced by the monomer inputs. |
| [monomer/02_TriP_S0.xyz](monomer/02_TriP_S0.xyz) | Optimized monomer geometry for the state-specific calculation. |
| [monomer/02_TriP_S1.xyz](monomer/02_TriP_S1.xyz) | Optimized monomer geometry for the state-specific calculation. |
| [monomer/03_MolX_D+@GS0_C1_desym.RasOrb](monomer/03_MolX_D+@GS0_C1_desym.RasOrb) | D+ monomer RasOrb at GS0. |
| [monomer/03_MolX_D+@GS0_C1_desym.in](monomer/03_MolX_D+@GS0_C1_desym.in) | D+ at GS0: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_D+@GS0_C1_desym.out](monomer/03_MolX_D+@GS0_C1_desym.out) | D+ at GS0: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_D+@GS1_C1_desym.RasOrb](monomer/03_MolX_D+@GS1_C1_desym.RasOrb) | D+ monomer RasOrb at GS1. |
| [monomer/03_MolX_D+@GS1_C1_desym.in](monomer/03_MolX_D+@GS1_C1_desym.in) | D+ at GS1: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_D+@GS1_C1_desym.out](monomer/03_MolX_D+@GS1_C1_desym.out) | D+ at GS1: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_D-@GS0_C1_desym.RasOrb](monomer/03_MolX_D-@GS0_C1_desym.RasOrb) | D- monomer RasOrb at GS0. |
| [monomer/03_MolX_D-@GS0_C1_desym.in](monomer/03_MolX_D-@GS0_C1_desym.in) | D- at GS0: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_D-@GS0_C1_desym.out](monomer/03_MolX_D-@GS0_C1_desym.out) | D- at GS0: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_D-@GS1_C1_desym.RasOrb](monomer/03_MolX_D-@GS1_C1_desym.RasOrb) | D- monomer RasOrb at GS1. |
| [monomer/03_MolX_D-@GS1_C1_desym.in](monomer/03_MolX_D-@GS1_C1_desym.in) | D- at GS1: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_D-@GS1_C1_desym.out](monomer/03_MolX_D-@GS1_C1_desym.out) | D- at GS1: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_S0@GS0_C1_desym.RasOrb](monomer/03_MolX_S0@GS0_C1_desym.RasOrb) | S0 monomer RasOrb at GS0. |
| [monomer/03_MolX_S0@GS0_C1_desym.in](monomer/03_MolX_S0@GS0_C1_desym.in) | S0 at GS0: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_S0@GS0_C1_desym.out](monomer/03_MolX_S0@GS0_C1_desym.out) | S0 at GS0: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_S0@GS1_C1_desym.RasOrb](monomer/03_MolX_S0@GS1_C1_desym.RasOrb) | S0 monomer RasOrb at GS1. |
| [monomer/03_MolX_S0@GS1_C1_desym.in](monomer/03_MolX_S0@GS1_C1_desym.in) | S0 at GS1: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_S0@GS1_C1_desym.out](monomer/03_MolX_S0@GS1_C1_desym.out) | S0 at GS1: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_S1@GS0_C1_desym.RasOrb](monomer/03_MolX_S1@GS0_C1_desym.RasOrb) | S1 monomer RasOrb at GS0. |
| [monomer/03_MolX_S1@GS0_C1_desym.in](monomer/03_MolX_S1@GS0_C1_desym.in) | S1 at GS0: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_S1@GS0_C1_desym.out](monomer/03_MolX_S1@GS0_C1_desym.out) | S1 at GS0: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_S1@GS1_C1_desym.RasOrb](monomer/03_MolX_S1@GS1_C1_desym.RasOrb) | S1 monomer RasOrb at GS1. |
| [monomer/03_MolX_S1@GS1_C1_desym.in](monomer/03_MolX_S1@GS1_C1_desym.in) | S1 at GS1: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_S1@GS1_C1_desym.out](monomer/03_MolX_S1@GS1_C1_desym.out) | S1 at GS1: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_T1@GS0_C1_desym.RasOrb](monomer/03_MolX_T1@GS0_C1_desym.RasOrb) | T1 monomer RasOrb at GS0. |
| [monomer/03_MolX_T1@GS0_C1_desym.in](monomer/03_MolX_T1@GS0_C1_desym.in) | T1 at GS0: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_T1@GS0_C1_desym.out](monomer/03_MolX_T1@GS0_C1_desym.out) | T1 at GS0: CASSCF/CASPT2 energies. |
| [monomer/03_MolX_T1@GS1_C1_desym.RasOrb](monomer/03_MolX_T1@GS1_C1_desym.RasOrb) | T1 monomer RasOrb at GS1. |
| [monomer/03_MolX_T1@GS1_C1_desym.in](monomer/03_MolX_T1@GS1_C1_desym.in) | T1 at GS1: state-specific CASSCF/CASPT2 input. |
| [monomer/03_MolX_T1@GS1_C1_desym.out](monomer/03_MolX_T1@GS1_C1_desym.out) | T1 at GS1: CASSCF/CASPT2 energies. |
| [monomer/Pentacene_monomer_S0_CAS66.rasscf.molden](monomer/Pentacene_monomer_S0_CAS66.rasscf.molden) | Monomer orbital source for the SI active-space figures. |
| [occupation_string_pairs.csv](occupation_string_pairs.csv) | 70-pair table for all three systems: 28 INDO pairs and 21 pairs for each pentacene aggregate; Fig. 5 and SI S9. |
| [trimer/B1_TriPABC_GS1S0S0.xyz](trimer/B1_TriPABC_GS1S0S0.xyz) | Aggregate geometry, in angstroms. |
| [trimer/B2_COMMONORB_TriPABC_GS1S0S0_CB1e-5](trimer/B2_COMMONORB_TriPABC_GS1S0S0_CB1e-5) | Aggregate CommonOrb basis for the trimer. |
| [trimer/B2_TriPA_D+@GS1.det](trimer/B2_TriPA_D+@GS1.det) | Fragment determinant expansion. |
| [trimer/B2_TriPA_D+@GS1.vec](trimer/B2_TriPA_D+@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPA_D-@GS1.det](trimer/B2_TriPA_D-@GS1.det) | Fragment determinant expansion. |
| [trimer/B2_TriPA_D-@GS1.vec](trimer/B2_TriPA_D-@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPA_S0@GS1.det](trimer/B2_TriPA_S0@GS1.det) | Fragment determinant expansion. |
| [trimer/B2_TriPA_S0@GS1.vec](trimer/B2_TriPA_S0@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPA_S1@GS1.det](trimer/B2_TriPA_S1@GS1.det) | Fragment determinant expansion. |
| [trimer/B2_TriPA_S1@GS1.vec](trimer/B2_TriPA_S1@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPA_T1@GS1.det](trimer/B2_TriPA_T1@GS1.det) | Fragment determinant expansion. |
| [trimer/B2_TriPA_T1@GS1.vec](trimer/B2_TriPA_T1@GS1.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPB_D+@GS0.det](trimer/B2_TriPB_D+@GS0.det) | Fragment determinant expansion. |
| [trimer/B2_TriPB_D+@GS0.vec](trimer/B2_TriPB_D+@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPB_D-@GS0.det](trimer/B2_TriPB_D-@GS0.det) | Fragment determinant expansion. |
| [trimer/B2_TriPB_D-@GS0.vec](trimer/B2_TriPB_D-@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPB_S0@GS0.det](trimer/B2_TriPB_S0@GS0.det) | Fragment determinant expansion. |
| [trimer/B2_TriPB_S0@GS0.vec](trimer/B2_TriPB_S0@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPB_T1@GS0.det](trimer/B2_TriPB_T1@GS0.det) | Fragment determinant expansion. |
| [trimer/B2_TriPB_T1@GS0.vec](trimer/B2_TriPB_T1@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPC_S0@GS0.det](trimer/B2_TriPC_S0@GS0.det) | Fragment determinant expansion. |
| [trimer/B2_TriPC_S0@GS0.vec](trimer/B2_TriPC_S0@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/B2_TriPC_T1@GS0.det](trimer/B2_TriPC_T1@GS0.det) | Fragment determinant expansion. |
| [trimer/B2_TriPC_T1@GS0.vec](trimer/B2_TriPC_T1@GS0.vec) | Fragment orbital coefficients in the aggregate CommonOrb basis. |
| [trimer/TriPABC.one](trimer/TriPABC.one) | One-electron integrals and nuclear repulsion. |
| [trimer/TriPABC_GS1S0S0_CB1e-5_CHO_Low_tau1e-5.json](trimer/TriPABC_GS1S0S0_CB1e-5_CHO_Low_tau1e-5.json) | Raw/normalized matrices, energies, CASPT2 corrections and timings. |
| [trimer/TriPABC_GS1S0S0_CB1e-5_CHO_Low_tau1e-5.occupation_statistics.json](trimer/TriPABC_GS1S0S0_CB1e-5_CHO_Low_tau1e-5.occupation_statistics.json) | 21-pair occupation-string statistics for Fig. 5 and SI S9. |
| [trimer/TriPABC_GS1S0S0_CB1e-5_CHO_Low_tau1e-5.sec](trimer/TriPABC_GS1S0S0_CB1e-5_CHO_Low_tau1e-5.sec) | Calculation event log. |
| [trimer/trimer.yaml](trimer/trimer.yaml) | NOCI input, with both screening thresholds equal to 1e-5. |
| [trimer/trimer_caspt2_shifted_hamiltonian_Eh.csv](trimer/trimer_caspt2_shifted_hamiltonian_Eh.csv) | CASPT2-shifted normalized Hamiltonian in Hartree; Fig. 6. |
| [trimer/trimer_normalized_hamiltonian_Eh.csv](trimer/trimer_normalized_hamiltonian_Eh.csv) | Unshifted normalized Hamiltonian in Hartree; SI S10. |
| [trimer/trimer_normalized_overlap.csv](trimer/trimer_normalized_overlap.csv) | Normalized overlap matrix; SI S10. |
