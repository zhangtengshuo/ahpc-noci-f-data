# Ethylene AHPC and reuse examples

SI S2/S3 and Fig. S1.

`reuse/` contains the complete CAS occupation-pattern enumeration for Fig. S1: 16 S1S0 determinants, six T1T1 determinants, 96 determinant pairs, 192 single-spin pair occurrences and 24 unique pairs.
The three CSV files give the determinant occupations, pair-to-string mapping and shared pair labels.
Regenerate them with `python reuse/enumerate_reuse.py OUTPUT_DIRECTORY`.
This enumeration uses the complete occupation patterns stated in the figure caption; it does not assign molecular CI coefficients.

`inputs/` contains the molecular NOCI inputs, fragment wave functions, geometry, CommonOrb, one-electron integrals, packed two-electron integrals, Cholesky factors, monomer RasOrb files and reference matrices.
The molecular calculation uses the fragment CI expansions in these input files.
`ethene_ahpc_complete.yaml` contains its AHPC intermediate quantities and determinant-level data.
The shared spatial overlap is $18\times18$ and its shared effective active-space overlap is $4\times4$.
`ethene_ahpc_dump.txt` records occupation-string mappings, counts and intermediate matrices.

Run the molecular calculation with `julia --project=/path/to/NOCI.jl ../run_noci.jl inputs/ahpc_cholesky.yaml OUTPUT_DIRECTORY ethene`.
Recalculate its reuse counts with `julia --project=/path/to/NOCI.jl ../reuse_counts.jl inputs/ahpc_cholesky.yaml OUTPUT.csv`.
Regenerate the intermediate quantities with `julia --project=/path/to/NOCI.jl ethene_ahpc_complete.jl OUTPUT_DIRECTORY`.
Run the detailed dump with `julia --project=/path/to/NOCI.jl ethene_ahpc_dump.jl inputs/ahpc_packed.yaml inputs/ahpc_cholesky.yaml inputs/conventional.matrix_snapshot.jld2 OUTPUT_DIRECTORY`.
Software setup is described in [SOFTWARE.md](../SOFTWARE.md).

## File list

| File | Contents / use |
|---|---|
| [B2_MolAB_CB.input](B2_MolAB_CB.input) | CommonOrb preparation input. |
| [B2_MolAB_CB.output](B2_MolAB_CB.output) | CommonOrb log confirming the retained 19+19=38 functions. |
| [ethene_ahpc_complete.jl](ethene_ahpc_complete.jl) | AHPC intermediate and matrix calculation script. |
| [ethene_ahpc_complete.yaml](ethene_ahpc_complete.yaml) | Readable AHPC data: orbital/overlap matrices and occupation-selected Schur blocks for SI S2. |
| [ethene_ahpc_dump.jl](ethene_ahpc_dump.jl) | Occupation-string, reuse and matrix output script. |
| [ethene_ahpc_dump.txt](ethene_ahpc_dump.txt) | Reference occupation-string mappings, reuse counts and AHPC intermediate quantities. |
| [inputs/03_MolX_D+@GS0_C1_desym.RasOrb](inputs/03_MolX_D+@GS0_C1_desym.RasOrb) | D+ monomer RasOrb. |
| [inputs/03_MolX_D-@GS0_C1_desym.RasOrb](inputs/03_MolX_D-@GS0_C1_desym.RasOrb) | D- monomer RasOrb. |
| [inputs/03_MolX_S0@GS0_C1_desym.RasOrb](inputs/03_MolX_S0@GS0_C1_desym.RasOrb) | S0 monomer RasOrb. |
| [inputs/03_MolX_S1@GS0_C1_desym.RasOrb](inputs/03_MolX_S1@GS0_C1_desym.RasOrb) | S1 monomer RasOrb. |
| [inputs/03_MolX_T1@GS0_C1_desym.RasOrb](inputs/03_MolX_T1@GS0_C1_desym.RasOrb) | T1 monomer RasOrb. |
| [inputs/B2_EtheneA_D+@GS0.det](inputs/B2_EtheneA_D+@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneA_D+@GS0.vec](inputs/B2_EtheneA_D+@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneA_D-@GS0.det](inputs/B2_EtheneA_D-@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneA_D-@GS0.vec](inputs/B2_EtheneA_D-@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneA_S0@GS0.det](inputs/B2_EtheneA_S0@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneA_S0@GS0.vec](inputs/B2_EtheneA_S0@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneA_S1@GS0.det](inputs/B2_EtheneA_S1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneA_S1@GS0.vec](inputs/B2_EtheneA_S1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneA_T1@GS0.det](inputs/B2_EtheneA_T1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneA_T1@GS0.vec](inputs/B2_EtheneA_T1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneB_D+@GS0.det](inputs/B2_EtheneB_D+@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneB_D+@GS0.vec](inputs/B2_EtheneB_D+@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneB_D-@GS0.det](inputs/B2_EtheneB_D-@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneB_D-@GS0.vec](inputs/B2_EtheneB_D-@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneB_S0@GS0.det](inputs/B2_EtheneB_S0@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneB_S0@GS0.vec](inputs/B2_EtheneB_S0@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneB_S1@GS0.det](inputs/B2_EtheneB_S1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneB_S1@GS0.vec](inputs/B2_EtheneB_S1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/B2_EtheneB_T1@GS0.det](inputs/B2_EtheneB_T1@GS0.det) | Fragment determinant expansion. |
| [inputs/B2_EtheneB_T1@GS0.vec](inputs/B2_EtheneB_T1@GS0.vec) | Fragment orbital coefficients in the dimer CommonOrb basis. |
| [inputs/COMMONORB](inputs/COMMONORB) | Ethylene dimer CommonOrb basis. |
| [inputs/EtheneAB.one](inputs/EtheneAB.one) | One-electron integrals and nuclear repulsion. |
| [inputs/EtheneAB.xyz](inputs/EtheneAB.xyz) | Ethylene dimer coordinates in angstroms. |
| [inputs/EtheneAB_001.two](inputs/EtheneAB_001.two) | Packed two-electron integrals. |
| [inputs/ahpc_cholesky.yaml](inputs/ahpc_cholesky.yaml) | NOCI input for the Cholesky calculation. |
| [inputs/ahpc_packed.yaml](inputs/ahpc_packed.yaml) | NOCI input for the packed-integral calculation. |
| [inputs/conventional.matrix_snapshot.jld2](inputs/conventional.matrix_snapshot.jld2) | Reference matrices read by the intermediate-output script. |
| [inputs/ethene_cho_factor_bundle.jld2](inputs/ethene_cho_factor_bundle.jld2) | Cholesky factor bundle. |
| [reuse/determinants.csv](reuse/determinants.csv) | The 16 S1S0 and six T1T1 determinant occupation patterns. |
| [reuse/enumerate_reuse.py](reuse/enumerate_reuse.py) | Enumerates the complete occupation-pattern reuse example of SI Fig. S1. |
| [reuse/pair_mapping.csv](reuse/pair_mapping.csv) | The 96 determinant-pair mappings to alpha/beta pair identifiers. |
| [reuse/unique_pairs.csv](reuse/unique_pairs.csv) | The 24 shared single-spin pair labels and their occurrence counts. |
