using NOCI
include(joinpath(@__DIR__, "dataset_input.jl"))

length(ARGS) in (2, 3) || error("Usage: run_noci.jl INPUT.yaml OUTPUT_DIRECTORY [LABEL]")
input = dataset_input(abspath(ARGS[1]))
label = length(ARGS) == 3 ? ARGS[3] : splitext(basename(ARGS[1]))[1]
NOCI.run_noci(input; output=NOCI.RunOutputOptions(abspath(ARGS[2]), label))
