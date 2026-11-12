using NOCI
using YAML
using LinearAlgebra

"""Resolve archived input references to files distributed beside each case."""
function dataset_file(config_path, value)
    name = basename(String(value))
    folder = dirname(abspath(config_path))
    candidates = [joinpath(folder, name), joinpath(folder, "inputs", name),
        joinpath(folder, "..", "inputs", name), joinpath(folder, "..", name)]
    for candidate in candidates
        isfile(candidate) && return abspath(candidate)
    end
    error("Required input is not distributed with this case: $name")
end

function dataset_pairwise_couplings(hamiltonian, overlap)
    n = size(hamiltonian, 1)
    couplings = zeros(eltype(hamiltonian), n, n)
    for i in 1:n, j in 1:i-1
        denominator = 1 - overlap[i,j]^2
        denominator > 0 || error("Linearly dependent state pair: $i, $j")
        value = (hamiltonian[i,j] - (hamiltonian[i,i] + hamiltonian[j,j]) * overlap[i,j] / 2) / denominator
        couplings[i,j] = couplings[j,i] = value
    end
    return couplings
end

function dataset_full_lowdin_couplings(hamiltonian, overlap)
    values, vectors = eigen(Symmetric(overlap))
    minimum(values) > 0 || error("The overlap matrix must be positive definite.")
    inverse_root = vectors * Diagonal(values .^ -0.5) * transpose(vectors)
    couplings = Matrix(Symmetric(transpose(inverse_root) * hamiltonian * inverse_root))
    for i in axes(couplings, 1)
        couplings[i,i] = 0
    end
    return couplings
end

function dataset_input(config_path; require_integrals=true)
    doc = YAML.load_file(config_path)
    for fragment in doc["system"]["fragment_states"]
        for key in ("determinants", "orbitals")
            fragment[key] = dataset_file(config_path, fragment[key])
        end
    end
    integrals = doc["integrals"]
    integrals["one_electron"] = dataset_file(config_path, integrals["one_electron"])
    if require_integrals
        if haskey(integrals, "two_electron")
            integrals["two_electron"] = [dataset_file(config_path, p) for p in integrals["two_electron"]]
        elseif haskey(integrals, "factor_bundle")
            integrals["factor_bundle"] = dataset_file(config_path, integrals["factor_bundle"])
        else
            error("This input requires two-electron data not distributed with the case.")
        end
    else
        doc["runtime"] = Dict("use_gpu" => false)
    end
    return mktempdir() do folder
        path = joinpath(folder, "input.yaml")
        YAML.write_file(path, doc)
        NOCI.read_noci_input(path)
    end
end
