# Pentacene dimer: software-release run examples

These files record the calculation inputs and outputs obtained with NOCI.jl 26.10.0.
The calculations use Julia 1.12.6 and OpenMolcas4NOCI 26.06.2.
Stages 01–04 contain the monomer states, aggregate geometry, CommonOrb basis and integral records used by the release-version NOCI calculations in stage 05.
The monomer wave functions use CAS(6,6) with the ANO-L-VDZP basis.
CASPT2 uses an IPEA shift of 0.0 Hartree and an imaginary shift of 0.2 Hartree.

| Workflow directory | Contents |
| --- | --- |
| `00_system_define/` | System definitions, geometries, initial orbitals and preparation inputs. |
| `01_fragment_states/` | Monomer CASSCF/CASPT2 inputs, outputs, orbitals and CI expansions. |
| `02_fragment_placement/` | Aggregate placement geometry. |
| `03_common_orbitals/` | CommonOrb records and transformed fragment wave functions. |
| `04_integrals/` | Integral records and conversion inputs/logs. |
| `05_noci/` | AHPC-NOCI inputs, Hamiltonian and overlap matrices, energies and timing records. |

| Run directory under `05_noci/` | Resources |
| --- | --- |
| `gpu8` | 8 physical CPU cores and 8 GPUs |

The reused CHO_Low integrals are RICD/acCD with actual SEWARD `CD Threshold: 1.00E-04`; CommonOrb and both NOCI screening thresholds are $10^{-5}$.
The dimer retains the GS1/GS0 geometry and includes the new conversion of the existing eight SEWARD shards to a factor bundle.

`results/result.out` and `results/result.json` contain the Hamiltonian and overlap matrices, NOCI energies and CASPT2 energy corrections.
`results/result.matrix_snapshot.jld2` stores the NOCI matrix snapshot; `run.time` and `results/result.sec` provide timing records.

The CommonOrb basis is supplied in `03_common_orbitals/COMMONORB`.
[Integral files not included](../OMITTED_FILES.md).
