# Benzene-dimer validation

Main-text Section 4.1; SI S6/S11.

## Software-release run examples

[release_version_run_examples/](release_version_run_examples/README.md) contains the inputs and outputs for calculations with NOCI.jl 26.10.0, including monomer states, CommonOrb bases and NOCI matrices and energies.
The existing `ahpc_packed/`, `ahpc_cholesky/`, `gronor/`, `inputs/` and `monomer/` contain the original datasets from earlier development versions; the remaining sections of this README describe those earlier data.

## Original development datasets

`gronor/Bz.out` and the two AHPC result-text files contain matrices, eigenvalues and determinant-pair counts.
The MEBF order is S0S0, S1S0, T1T1, D+D-, D-D+.
Hamiltonian elements and energies are in Hartree; overlaps are dimensionless.
The prefixes `raw_*` and `norm_*` identify unnormalized and normalized matrices.

`inputs/` contains fragment wave functions, geometry, CommonOrb and one-electron integrals.
`monomer/` contains the five state-specific RasOrb files and the inputs/outputs at the S0 geometry.
The monomer inputs use D2h symmetry followed by desymmetrization.
The packed two-electron file and Cholesky factor bundle are not included.

Analyze an output with `python ../analyze_matrices.py ahpc_packed/ahpc_oldchain_result.txt OUTPUT_DIRECTORY`.
Recalculate reuse counts with `julia --project=/path/to/NOCI.jl ../reuse_counts.jl ahpc_packed/ahpc_oldchain.yaml OUTPUT.csv`.
Software setup is described in [SOFTWARE.md](../SOFTWARE.md).

## File list

| File | Contents / use |
|---|---|
| [ahpc_cholesky/ahpc_cho_gs0s0.yaml](ahpc_cholesky/ahpc_cho_gs0s0.yaml) | Calculation input. |
| [ahpc_cholesky/ahpc_cho_result.txt](ahpc_cholesky/ahpc_cho_result.txt) | Matrix/eigenvalue/count output. |
| [ahpc_packed/ahpc_oldchain.yaml](ahpc_packed/ahpc_oldchain.yaml) | Calculation input. |
| [ahpc_packed/ahpc_oldchain_result.txt](ahpc_packed/ahpc_oldchain_result.txt) | Matrix/eigenvalue/count output. |
| [ahpc_packed/ahpc_run.jl](ahpc_packed/ahpc_run.jl) | Matrix and determinant-pair output driver. |
| [gronor/Bz.inp](gronor/Bz.inp) | GronOR benchmark input. |
| [gronor/Bz.out](gronor/Bz.out) | GronOR output: matrices, eigenvalues and nullity counts. |
| [inputs/B1_BzAB_GS0S0.xyz](inputs/B1_BzAB_GS0S0.xyz) | S0/S0 dimer geometry, in angstroms. |
| [inputs/B2_COMMONORB_BzAB_GS0S0_CB1e-4](inputs/B2_COMMONORB_BzAB_GS0S0_CB1e-4) | CommonOrb basis for the benzene dimer. |
| [inputs/BzAB.one](inputs/BzAB.one) | One-electron integrals and nuclear repulsion. |
| [inputs/BzAB.sys](inputs/BzAB.sys) | Native atom/basis system record. |
| [inputs/BzA_D+.det](inputs/BzA_D+.det) | Fragment determinant expansion. |
| [inputs/BzA_D+.vec](inputs/BzA_D+.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzA_D-.det](inputs/BzA_D-.det) | Fragment determinant expansion. |
| [inputs/BzA_D-.vec](inputs/BzA_D-.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzA_S0.det](inputs/BzA_S0.det) | Fragment determinant expansion. |
| [inputs/BzA_S0.vec](inputs/BzA_S0.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzA_S1.det](inputs/BzA_S1.det) | Fragment determinant expansion. |
| [inputs/BzA_S1.vec](inputs/BzA_S1.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzA_T1.det](inputs/BzA_T1.det) | Fragment determinant expansion. |
| [inputs/BzA_T1.vec](inputs/BzA_T1.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzB_D+.det](inputs/BzB_D+.det) | Fragment determinant expansion. |
| [inputs/BzB_D+.vec](inputs/BzB_D+.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzB_D-.det](inputs/BzB_D-.det) | Fragment determinant expansion. |
| [inputs/BzB_D-.vec](inputs/BzB_D-.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzB_S0.det](inputs/BzB_S0.det) | Fragment determinant expansion. |
| [inputs/BzB_S0.vec](inputs/BzB_S0.vec) | Fragment CommonOrb orbital coefficients. |
| [inputs/BzB_T1.det](inputs/BzB_T1.det) | Fragment determinant expansion. |
| [inputs/BzB_T1.vec](inputs/BzB_T1.vec) | Fragment CommonOrb orbital coefficients. |
| [monomer/01_MolX_CAS66_D2h.RasOrb](monomer/01_MolX_CAS66_D2h.RasOrb) | Initial RasOrb used by the five monomer inputs. |
| [monomer/02CAS_Bz_S0_opt.Opt.xyz](monomer/02CAS_Bz_S0_opt.Opt.xyz) | Optimized S0 monomer coordinates referenced by all five fragment inputs. |
| [monomer/03_MolX_D+@GS0_Ci_desym.RasOrb](monomer/03_MolX_D+@GS0_Ci_desym.RasOrb) | D+ monomer RasOrb at the S0 geometry. |
| [monomer/03_MolX_D+@GS0_Ci_desym.in](monomer/03_MolX_D+@GS0_Ci_desym.in) | D+ fragment OpenMolcas input at the ground-state geometry. |
| [monomer/03_MolX_D+@GS0_Ci_desym.out](monomer/03_MolX_D+@GS0_Ci_desym.out) | D+ fragment OpenMolcas output. |
| [monomer/03_MolX_D-@GS0_Ci_desym.RasOrb](monomer/03_MolX_D-@GS0_Ci_desym.RasOrb) | D- monomer RasOrb at the S0 geometry. |
| [monomer/03_MolX_D-@GS0_Ci_desym.in](monomer/03_MolX_D-@GS0_Ci_desym.in) | D- fragment OpenMolcas input at the ground-state geometry. |
| [monomer/03_MolX_D-@GS0_Ci_desym.out](monomer/03_MolX_D-@GS0_Ci_desym.out) | D- fragment OpenMolcas output. |
| [monomer/03_MolX_S0@GS0_Ci_desym.RasOrb](monomer/03_MolX_S0@GS0_Ci_desym.RasOrb) | S0 monomer RasOrb at the S0 geometry. |
| [monomer/03_MolX_S0@GS0_Ci_desym.in](monomer/03_MolX_S0@GS0_Ci_desym.in) | S0 fragment OpenMolcas input at the ground-state geometry. |
| [monomer/03_MolX_S0@GS0_Ci_desym.out](monomer/03_MolX_S0@GS0_Ci_desym.out) | S0 fragment OpenMolcas output. |
| [monomer/03_MolX_S0@GS0_Ci_desym.rasscf.molden](monomer/03_MolX_S0@GS0_Ci_desym.rasscf.molden) | S0 active-space orbitals for the SI figure. |
| [monomer/03_MolX_S1@GS0_Ci_desym.RasOrb](monomer/03_MolX_S1@GS0_Ci_desym.RasOrb) | S1 monomer RasOrb at the S0 geometry. |
| [monomer/03_MolX_S1@GS0_Ci_desym.in](monomer/03_MolX_S1@GS0_Ci_desym.in) | S1 fragment OpenMolcas input at the ground-state geometry. |
| [monomer/03_MolX_S1@GS0_Ci_desym.out](monomer/03_MolX_S1@GS0_Ci_desym.out) | S1 fragment OpenMolcas output. |
| [monomer/03_MolX_T1@GS0_Ci_desym.RasOrb](monomer/03_MolX_T1@GS0_Ci_desym.RasOrb) | T1 monomer RasOrb at the S0 geometry. |
| [monomer/03_MolX_T1@GS0_Ci_desym.in](monomer/03_MolX_T1@GS0_Ci_desym.in) | T1 fragment OpenMolcas input at the ground-state geometry. |
| [monomer/03_MolX_T1@GS0_Ci_desym.out](monomer/03_MolX_T1@GS0_Ci_desym.out) | T1 fragment OpenMolcas output. |
