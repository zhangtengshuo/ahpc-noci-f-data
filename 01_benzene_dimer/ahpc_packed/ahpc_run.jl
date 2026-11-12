using NOCI, Printf

const yaml_path = abspath(ARGS[1])
const out_path = abspath(ARGS[2])

input = read_noci_input(yaml_path)
result = evaluate(input)

function dump_matrix(io, name, M)
    for i in axes(M, 1)
        println(io, name, "|", join((@sprintf "%23.15e" M[i, j]) for j in axes(M, 2)), "|")
    end
end

open(out_path, "w") do io
    r = result.raw
    println(io, "labels|", join(r.labels, ","))
    dump_matrix(io, "raw_overlap", r.overlap)
    dump_matrix(io, "raw_h1_electronic", r.electronic_one_body)
    dump_matrix(io, "raw_nuclear", r.nuclear_repulsion)
    dump_matrix(io, "raw_h1_with_nuclear", r.one_body_with_nuclear_repulsion)
    dump_matrix(io, "raw_h2", r.two_body)
    dump_matrix(io, "raw_hamiltonian", r.hamiltonian)
    dump_matrix(io, "norm_overlap", result.normalized.overlap)
    dump_matrix(io, "norm_h1_with_nuclear", result.normalized.one_body_with_nuclear_repulsion)
    dump_matrix(io, "norm_h2", result.normalized.two_body)
    dump_matrix(io, "norm_hamiltonian", result.normalized.hamiltonian)
    for i in eachindex(r.labels)
        println(io, "paircount|", r.labels[i], "|", join(r.determinant_pair_counts[i, :], ","))
    end
    println(io, "total_pairs|", sum(r.determinant_pair_counts))
    println(io, "energies|", join((@sprintf "%.15e" e for e in result.solution.energies), ","))
    for (k, v) in sort(collect(result.telemetry.counters))
        println(io, "counter|", k, "|", v)
    end
end
println("AHPC_DONE total_pairs=", sum(result.raw.determinant_pair_counts))
