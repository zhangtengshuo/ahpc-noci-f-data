#!/usr/bin/env julia
# Calculate the AHPC intermediate quantities for the ethylene example.
# Usage: julia --project=/path/to/NOCI.jl ethene_ahpc_complete.jl OUTPUT_DIRECTORY [INPUT.yaml]
using NOCI, LinearAlgebra, JLD2, YAML

length(ARGS) in (1, 2) || error("Usage: ethene_ahpc_complete.jl OUTPUT_DIRECTORY [INPUT.yaml]")
outdir = abspath(ARGS[1])
mkpath(outdir)
include(joinpath(@__DIR__, "..", "dataset_input.jl"))
config = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(@__DIR__, "inputs", "ahpc_cholesky.yaml")
input = dataset_input(config)
problem = NOCI.load_problem_definition(input)
bfs = NOCI.build_molecular_basis_functions(problem; screening=input.screening)
one = NOCI.read_one_electron_file(input.backend.one_path).values
factors = NOCI.read_cholesky_factor_bundle(input.backend.source.factor_path).factors
runtime = NOCI.RuntimeOptions()
canonical = NOCI.build_canonical_mebf_spin_strings.(bfs)
mat = NOCI.assemble_ahpc_spin_string_matrices(bfs, one, factors, input.screening)
result = NOCI.evaluate(input)
shifts = [b.correlation_energy_shift for b in bfs]
shifted_hamiltonian = NOCI.apply_diagonal_shifts(result.normalized.hamiltonian, result.normalized.overlap, shifts; weighting=:identity)
shifted_eigen = eigen(Symmetric(shifted_hamiltonian), Symmetric(result.normalized.overlap))

"""Convert numerical arrays to portable row-oriented YAML without reducing precision."""
plain(x::AbstractMatrix) = [collect(x[i, :]) for i in axes(x, 1)]
plain(x::AbstractVector) = [plain(v) for v in x]
plain(x::NamedTuple) = Dict(string(k) => plain(v) for (k, v) in pairs(x))
plain(x::AbstractDict) = Dict(string(k) => plain(v) for (k, v) in x)
plain(x) = x

"""Return occupied active-pool positions from a canonical spin string."""
positions(bits) = NOCI.spin_string_indices(bits, 4)

"""Evaluate the n=2 PDDC outputs and check the production polynomial routine."""
function local_objects(prep, coefficients, leftbits, rightbits)
    t, u = positions(leftbits), positions(rightbits)
    X = prep.effective_active_overlap[t, u]
    a, c, b, d = vec(X)
    C0 = a*d-b*c
    C1 = [d -c; -b a]
    xi = [t[i]+4*(u[j]-1) for j in 1:2 for i in 1:2]
    R = (coefficients.coulomb_response-coefficients.exchange_response)[xi, xi]
    R2 = R[1,4]+R[4,1]-R[3,2]-R[2,3]
    polynomial = NOCI.polynomial_contracted_complementary_cofactors(X, R)
    @assert isapprox(polynomial.cofactor_scalar, C0; atol=2e-15)
    @assert maximum(abs.(polynomial.first_cofactor-C1)) < 2e-15
    @assert isapprox(polynomial.contracted_second, R2; atol=2e-14)
    return (t=t, u=u, xi=xi, X=X, C0=C0, C1=C1, R=R, R2=R2,
        l1=dot(coefficients.one_body_active[xi], vec(C1)),
        l2=dot(coefficients.two_body_inactive_active[xi], vec(C1)))
end

blocks = Dict{String, Any}()
maxerror = zeros(4)
for li in eachindex(bfs), ri in 1:li
    L, R = bfs[li], bfs[ri]
    prep = NOCI.prepare_ahpc_state_pair(L, R, one.overlap; tau_singularity=input.screening.tau_singularity)
    coeff = NOCI.precontract_active_hamiltonian(one.kinetic+one.nuclear_attraction, prep.response_space, factors; runtime)
    pairlist = NOCI._build_ahpc_state_pair_pairs(L, R, input.screening, li, ri; runtime, telemetry=nothing)
    rows = Any[]
    totals = zeros(4)
    for I in 1:pairlist.n_left, J in (pairlist.diagonal ? I : 1):pairlist.last_right[I]
        A = local_objects(prep, coeff, L.expansion.active_occupations.alpha_bits[I], R.expansion.active_occupations.alpha_bits[J])
        B = local_objects(prep, coeff, L.expansion.active_occupations.beta_bits[I], R.expansion.active_occupations.beta_bits[J])
        q = prep.inactive_overlap_determinant^2
        samealpha = 0.25*A.R2*B.C0
        samebeta = 0.25*A.C0*B.R2
        opposite = dot(vec(A.C1), coeff.coulomb_response[A.xi, B.xi]*vec(B.C1))
        terms1 = [2*coeff.one_body_inactive*A.C0*B.C0, A.l1*B.C0, A.C0*B.l1]
        terms2 = [coeff.two_body_inactive_inactive*A.C0*B.C0, A.l2*B.C0, A.C0*B.l2, samealpha, samebeta, opposite]
        S = q*A.C0*B.C0
        H1 = q*sum(terms1)
        H2 = q*sum(terms2)
        Hnn = one.nuclear_repulsion*S
        phase = canonical[li].determinant_phase[I]*canonical[ri].determinant_phase[J]
        multiplicity = li==ri && I!=J ? 2 : 1
        weight = multiplicity*L.expansion.coefficients[I]*R.expansion.coefficients[J]*phase
        vals = [S, H1+Hnn, H2, H1+Hnn+H2]
        totals .+= weight.*vals
        push!(rows, (I=I, J=J, weight=weight, phase=phase, multiplicity=multiplicity,
            alpha=A, beta=B, S=S, H1=H1, Hnn=Hnn, H2=H2, H=sum(vals[2:3]),
            terms1=terms1, terms2=terms2, weighted=weight.*vals))
    end
    target = [mat.overlap[li,ri], mat.one_body_with_nuclear_repulsion[li,ri], mat.two_body[li,ri], mat.hamiltonian[li,ri]]
    maxerror .= max.(maxerror, abs.(totals-target))
    maximum(abs.(totals-target)) > 5e-10 && println("Mismatch ", li, " ", ri, ": ", totals-target)
    blocks["$(li)_$(ri)"] = (left=L.label, right=R.label, preparation=prep,
        coefficients=coeff, rows=rows, totals=totals, target=target)
end
println("Assembly check: ", maxerror)
@assert maxerror[1]<1e-12 && maximum(maxerror[2:4])<5e-10

# Portable full-precision supplement: all 253 retained pairs and their PDDC data.
tt = blocks["4_2"]
prep, coeff = tt.preparation, tt.coefficients
space = prep.response_space
J_II, K_II = 0.0, 0.0
J_IA, K_IA = zeros(16), zeros(16)
projection_example = nothing
factor_matrix = zeros(38,38)
for auxiliary in axes(factors.coefficients, 2)
    NOCI._unpack_cholesky_factor!(factor_matrix, factors, auxiliary)
    zI = sum(factor_matrix.*space.inactive_density_kernel)
    zA = transpose(space.left_active_response)*factor_matrix*space.right_active_response
    kI = sum(space.inactive_density_kernel.*(factor_matrix*transpose(space.inactive_density_kernel)*factor_matrix))
    kA = transpose(space.left_active_response)*factor_matrix*transpose(space.inactive_density_kernel)*factor_matrix*space.right_active_response
    global J_II += zI*zI
    global K_II += kI
    J_IA .+= zI.*vec(zA)
    K_IA .+= vec(kA)
    auxiliary==1 && (global projection_example=(zI=zI, zA=zA, kI=kI, kA=kA))
end
@assert abs(2J_II-K_II-coeff.two_body_inactive_inactive)<1e-11
@assert maximum(abs.(2J_IA-K_IA-coeff.two_body_inactive_active))<1e-11
supplement = Dict(
    "configuration"=>config, "nuclear_repulsion"=>one.nuclear_repulsion,
    "max_assembly_errors"=>maxerror,
    "mebfs"=>[(label=b.label, coefficients=b.expansion.coefficients,
        orbital_coefficients=b.orbital_coefficients,
        alpha=[positions(v) for v in b.expansion.active_occupations.alpha_bits],
        beta=[positions(v) for v in b.expansion.active_occupations.beta_bits],
        phase=canonical[i].determinant_phase, shift=b.correlation_energy_shift) for (i,b) in enumerate(bfs)],
    "example"=>(det_II=prep.inactive_overlap_determinant, rcond=prep.inactive_reciprocal_condition,
        X_II=prep.X_inactive_inactive, X_IA=prep.X_inactive_active, X_AI=prep.X_active_inactive,
        X_AA=prep.X_active_active, Xtilde=prep.effective_active_overlap,
        C_L=bfs[4].orbital_coefficients, C_R=bfs[2].orbital_coefficients,
        Ctilde_L=prep.response_space.left_active_response, Ctilde_R=prep.response_space.right_active_response,
        P_I=prep.response_space.inactive_density_kernel,
        o_I=coeff.one_body_inactive, o_A=coeff.one_body_active,
        J_II=J_II, K_II=K_II, J_IA=J_IA, K_IA=K_IA,
        projection_first=projection_example,
        h_II=coeff.two_body_inactive_inactive, h_IA=coeff.two_body_inactive_active,
        J_AA=coeff.coulomb_response, K_AA=coeff.exchange_response,
        rows=tt.rows, totals=tt.totals),
    "blocks"=>Dict(k=>(left=v.left, right=v.right, rows=v.rows, totals=v.totals,
        Xtilde_AA=v.preparation.effective_active_overlap,
        det_II=v.preparation.inactive_overlap_determinant,
        rcond=v.preparation.inactive_reciprocal_condition) for (k,v) in blocks),
    "raw"=>(S=mat.overlap, H1nn=mat.one_body_with_nuclear_repulsion, H2=mat.two_body, H=mat.hamiltonian),
    "normalized"=>(S=result.normalized.overlap, H=result.normalized.hamiltonian),
    "energies"=>result.solution.energies, "vectors"=>result.solution.vectors,
    "couplings_pairwise"=>dataset_pairwise_couplings(result.normalized.hamiltonian, result.normalized.overlap),
    "couplings_full"=>dataset_full_lowdin_couplings(result.normalized.hamiltonian, result.normalized.overlap),
    "caspt2_identity"=>(shifts=shifts, H=shifted_hamiltonian, energies=shifted_eigen.values),
    "basis_overlap"=>one.overlap,
    "one_body_operator"=>one.kinetic+one.nuclear_attraction)
JLD2.save(joinpath(outdir, "ethene_ahpc_complete.jld2"), supplement)
YAML.write_file(joinpath(outdir, "ethene_ahpc_complete.yaml"), plain(supplement))
println("Verified all 21 blocks and 253 determinant pairs; maximum errors (S,H1+nn,H2,H): ", maxerror)
