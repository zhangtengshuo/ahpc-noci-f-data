using NOCI
include(joinpath(@__DIR__, "dataset_input.jl"))

length(ARGS) == 2 || error("Usage: reuse_counts.jl INPUT.yaml OUTPUT.csv")
input = dataset_input(abspath(ARGS[1]); require_integrals=false)
runtime = NOCI.RuntimeOptions()
problem = NOCI.load_problem_definition(input; runtime)
basis = NOCI.build_molecular_basis_functions(problem; screening=input.screening, runtime)
canonical = NOCI.build_canonical_mebf_spin_strings.(basis)
csv(value) = "\"" * replace(string(value), "\"" => "\"\"") * "\""
open(ARGS[2], "w") do io
    println(io, "state_left,state_right,determinant_pairs,alpha_pairs,beta_pairs,shared_pairs,union_pairs,reuse_factor")
    for i in eachindex(basis), j in 1:i
        block = NOCI._build_ahpc_state_pair_pairs(basis[i], basis[j], input.screening, i, j;
            runtime, telemetry=nothing)
        incidence = NOCI.build_spin_pair_incidence(basis[i], basis[j], canonical[i], canonical[j],
            block; support_lookup=NOCI.AdaptiveLookup, runtime, telemetry=nothing)
        support = incidence.support
        n = block.total_count
        u = length(support.union_support)
        values = (basis[i].label, basis[j].label, n, length(support.alpha_support),
            length(support.beta_support), support.intersection_count, u, u == 0 ? 0.0 : 2n/u)
        println(io, join(csv.(values), ","))
    end
end
