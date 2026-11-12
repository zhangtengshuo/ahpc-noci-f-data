# Pentacene: software-release run examples

This directory contains inputs and outputs for AHPC-NOCI calculations on pentacene dimers and trimers using NOCI.jl 26.10.0.
The monomer wave functions use CAS(6,6) with the ANO-L-VDZP basis.

- [dimer/](dimer/README.md): GS1/GS0 dimer, 8 physical CPU cores and 8 GPUs, including conversion of the existing SEWARD data.
- [trimer/](trimer/README.md): GS1/GS0/GS0 trimer, 8 physical CPU cores and 8 GPUs, reusing the existing factor bundle.

Each calculation includes monomer CASSCF/CASPT2 states, aggregate geometry, the CommonOrb basis, integral records and NOCI matrices.
The actual RICD/acCD integral threshold is $10^{-4}$, and the CommonOrb and NOCI screening thresholds are $10^{-5}$.


CommonOrb basis sets are supplied in both subdirectories.
[Integral files not included](OMITTED_FILES.md).
