#!/usr/bin/env julia
# =============================================================================
# ethene_ahpc_dump.jl — AHPC-NOCI intermediate-quantity dump for the ethene
# dimer demonstration (manuscript: AHPC-NOCI, spin_string branch).
#
# Reproduces the production code path
#   read_noci_input -> load_problem_definition -> build_molecular_basis_functions
#   -> prepare_ahpc_state_pair (X blocks, Schur-reduced X~_AA)
#   -> _build_ahpc_state_pair_pairs (screened determinant-pair block)
#   -> build_spin_pair_incidence (occupation-string-pair support U_alpha/U_beta/U)
#   -> assemble_ahpc_spin_string_matrices (raw NOCI matrices + telemetry)
#   -> evaluate (solution) + Loewdin couplings
# and prints every object in human-readable form, saving the same content to
# a JLD2 archive.
#
# usage:
#   julia --project=<path-to-NOCI.jl> ethene_ahpc_dump.jl \
#       <packed_yaml> <cholesky_yaml> <reference_jld2> <outdir>
# =============================================================================

using NOCI
using JLD2
using LinearAlgebra
using Printf
using Dates
include(joinpath(@__DIR__, "..", "dataset_input.jl"))

const PACKED_YAML = abspath(ARGS[1])
const CHO_YAML = abspath(ARGS[2])
const REF_JLD2 = abspath(ARGS[3])
const OUTDIR = abspath(ARGS[4])

mkpath(OUTDIR)
const TXT = joinpath(OUTDIR, "ethene_ahpc_dump.txt")
const io = open(TXT, "w")

pr(args...) = (s = string(args...); println(io, s); flush(io))

function prmat(title::String, M::AbstractMatrix; digits::Int=10, width::Int=18)
    pr(title, "  (", size(M, 1), "x", size(M, 2), ")")
    for i in axes(M, 1)
        row = "  "
        for j in axes(M, 2)
            row *= @sprintf(" %*.*f", width, digits, M[i, j])
        end
        pr(row)
    end
end

prvec(title::String, v::AbstractVector; digits::Int=10, width::Int=18) =
    pr(title, " ", join((@sprintf("%*.*f", width, digits, x) for x in v), ""))

# decode one active-orbital occupation digit: 2 double, a alpha, b beta, 0 empty
function occ_digit(alpha_bits::UInt64, beta_bits::UInt64, p::Int)
    a = (alpha_bits >> (p - 1)) & one(UInt64)
    b = (beta_bits >> (p - 1)) & one(UInt64)
    a == 1 && b == 1 && return '2'
    a == 1 && return 'a'
    b == 1 && return 'b'
    return '0'
end

bits_string(bits::UInt64, n::Int) = join((bits >> (p - 1)) & one(UInt64) for p in 1:n)

# =============================================================================
pr("# =========================================================================")
pr("# AHPC-NOCI intermediate dump — ethene dimer (CAS(2,2) x 2, cc-pVDZ, C1)")
pr("# generated: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
pr("# packed yaml : ", PACKED_YAML)
pr("# cholesky yaml: ", CHO_YAML)
pr("# =========================================================================")

# ---------------------------------------------------------------- 1. problem
input = dataset_input(PACKED_YAML)
screening = input.screening
runtime = NOCI.RuntimeOptions()
problem = NOCI.load_problem_definition(input)
bfs = NOCI.build_molecular_basis_functions(problem; screening=screening)
one_el = NOCI.read_one_electron_file(input.backend.one_path).values
cho_input = dataset_input(CHO_YAML)
factors = NOCI.read_cholesky_factor_bundle(cho_input.backend.source.factor_path).factors

n_mebf = length(bfs)
n_basis = size(one_el.overlap, 1)
pr("")
pr("# ---------------------------------------------------------------")
pr("# 1. problem dimensions")
pr("# ---------------------------------------------------------------")
pr("CommonOrb working basis dim (one-electron overlap): ", n_basis)
pr("basis-overlap eigenvalue range: ",
   minimum(eigvals(Symmetric(one_el.overlap))), " .. ", maximum(eigvals(Symmetric(one_el.overlap))))
pr("screening: tau_ci=", screening.tau_ci, "  tau_ci_offdiag=", screening.tau_ci_offdiagonal,
   "  tau_singularity=", screening.tau_singularity)
pr("Cholesky factor shape: ", size(factors.coefficients), " (canonical pairs x n_chol)")

# ------------------------------------------------------------- 2. MEBF table
pr("")
pr("# ---------------------------------------------------------------")
pr("# 2. MEBF determinant expansions (Eq. 1-3 level)")
pr("# ---------------------------------------------------------------")
canonical = NOCI.build_canonical_mebf_spin_strings.(bfs)
for (mu, b) in enumerate(bfs)
    ex = b.expansion
    n_act = size(ex.active_occupations, 1)
    pr("")
    pr("MEBF ", mu, ": ", b.label, "   fragment-state indices = ", b.fragment_state_indices)
    pr("  orbital_coefficients: ", size(b.orbital_coefficients), "  inactive = ", ex.inactive_orbitals,
       "  active = ", n_act, "  determinants = ", length(ex.coefficients))
    pr("  correlation_energy_shift = ", b.correlation_energy_shift)
    for d in eachindex(ex.coefficients)
        ab = ex.active_occupations.alpha_bits[d]
        bb = ex.active_occupations.beta_bits[d]
        perorb = join(occ_digit(ab, bb, p) for p in 1:n_act)
        pr(@sprintf("   det %2d  c = %+.10f  phase = %+d  occ = %s   alpha = (%s)  beta = (%s)",
                    d, ex.coefficients[d], canonical[mu].determinant_phase[d], perorb,
                    bits_string(ab, n_act), bits_string(bb, n_act)))
    end
end

# --------------------------------------------------- 3. per-MEBF-pair objects
pr("")
pr("# ---------------------------------------------------------------")
pr("# 3. state-pair preparations and occupation-string-pair support")
pr("# ---------------------------------------------------------------")

store = Dict{String,Any}()

for left_mebf in 1:n_mebf, right_mebf in 1:left_mebf
    L, R = bfs[left_mebf], bfs[right_mebf]
    pr("")
    pr("=========== MEBF pair (", L.label, ", ", R.label, ")  [", left_mebf, ",", right_mebf, "] ===========")
    pairs_block = NOCI._build_ahpc_state_pair_pairs(L, R, screening, left_mebf, right_mebf;
        runtime, telemetry=nothing)
    pr("retained determinant pairs (row-compressed): ", pairs_block.total_count,
       "   [n_left=", pairs_block.n_left, " n_right=", pairs_block.n_right,
       " diagonal=", pairs_block.diagonal, "]")

    prep = NOCI.prepare_ahpc_state_pair(L, R, one_el.overlap;
        tau_singularity=screening.tau_singularity, runtime, telemetry=nothing)
    pr("X_II det = ", prep.inactive_overlap_determinant,
       "   rcond = ", prep.inactive_reciprocal_condition)
    prmat("X_II (inactive x inactive):", prep.X_inactive_inactive; digits=8)
    prmat("X_AA (active pool, left rows / right cols):", prep.X_active_active; digits=8)
    prmat("X_AI:", prep.X_active_inactive; digits=8)
    prmat("X~_AA = X_AA - X_AI X_II^-1 X_IA (Schur):", prep.effective_active_overlap; digits=8)
    prmat("left active-response C_A^L (first 8 rows):", prep.response_space.left_active_response[1:min(8, end), :]; digits=6)

    incidence = NOCI.build_spin_pair_incidence(L, R, canonical[left_mebf], canonical[right_mebf],
        pairs_block; support_lookup=NOCI.AdaptiveLookup, runtime, telemetry=nothing)
    support = incidence.support
    n_act = size(L.expansion.active_occupations, 1)
    pr("support counts: |U_alpha| = ", length(support.alpha_support),
       "  |U_beta| = ", length(support.beta_support),
       "  |U| = ", length(support.union_support),
       "  intersection = ", support.intersection_count)
    pr("")
    pr("  union string pairs (rho): t = left occupied list, u = right occupied list,")
    pr("  X~_sigma = X~_AA[t,u], det, rank, omega = n_occ - rank (tau_sing = ", screening.tau_singularity, "):")
    for (uid, pair) in enumerate(support.union_support)
        t = NOCI.spin_string_indices(pair.left_bits, n_act)
        u = NOCI.spin_string_indices(pair.right_bits, n_act)
        Xs = prep.effective_active_overlap[t, u]
        sv = svdvals(Xs)
        rk = count(>(screening.tau_singularity * max(1.0, sv[1])), sv)
        omega = length(t) - rk
        pr(@sprintf("   rho %3d  alphaL=(%s) alphaR=(%s)  t=%s u=%s  n_occ=%d  det=%+ .3e  rank=%d  omega=%d  svmin=% .3e",
                    uid, bits_string(pair.left_bits, n_act), bits_string(pair.right_bits, n_act),
                    join(t, ","), join(u, ","), length(t), det(Xs), rk, omega,
                    isempty(sv) ? 0.0 : sv[end]))
        prmat("     X~_sigma:", Xs; digits=6, width=14)
    end

    store["pair_$(left_mebf)_$(right_mebf)"] = (
        label_left=L.label, label_right=R.label,
        n_pairs=pairs_block.total_count,
        X_II=prep.X_inactive_inactive, X_IA=prep.X_inactive_active,
        X_AI=prep.X_active_inactive, X_AA=prep.X_active_active,
        det_II=prep.inactive_overlap_determinant, rcond=prep.inactive_reciprocal_condition,
        Xtilde_AA=prep.effective_active_overlap,
        union_alpha=[(p.left_bits, p.right_bits) for p in support.alpha_support],
        union_beta=[(p.left_bits, p.right_bits) for p in support.beta_support],
        union_all=[(p.left_bits, p.right_bits) for p in support.union_support],
        intersection=support.intersection_count,
    )
end

# ------------------------------------ 3b. determinant-pair / weight demo for
#                                     the S1S0-T1T1 block 
pr("")
pr("# ---------------------------------------------------------------")
pr("# 3b. determinant-pair -> string-pair incidence (S1S0 / T1T1 block)")
pr("# ---------------------------------------------------------------")
li = findfirst(b -> b.label == "S1S0", bfs)
tj = findfirst(b -> b.label == "T1T1", bfs)
left_idx, right_idx = max(li, tj), min(li, tj)   # production enumerates lower-triangular blocks
L, R = bfs[left_idx], bfs[right_idx]
pr("block orientation: left = ", L.label, " (", left_idx, "), right = ", R.label, " (", right_idx, ")")
blk = NOCI._build_ahpc_state_pair_pairs(L, R, screening, left_idx, right_idx; runtime, telemetry=nothing)
n_act = size(L.expansion.active_occupations, 1)
occ = Dict{Tuple{UInt64,UInt64},Vector{Tuple{Int,Int,Float64}}}()
for I in 1:blk.n_left, J in 1:blk.last_right[I]
    akey = (L.expansion.active_occupations.alpha_bits[I], R.expansion.active_occupations.alpha_bits[J])
    bkey = (L.expansion.active_occupations.beta_bits[I], R.expansion.active_occupations.beta_bits[J])
    for key in (akey, bkey)
        push!(get!(occ, key, Tuple{Int,Int,Float64}[]), (I, J, L.expansion.coefficients[I] * R.expansion.coefficients[J]))
    end
end
pr("retained pairs = ", blk.total_count, "; spin-resolved occurrences = ", 2 * blk.total_count,
   "; distinct single-spin string pairs = ", length(occ))
multiplicities = Dict{Int,Int}()
for v in values(occ)
    multiplicities[length(v)] = get(multiplicities, length(v), 0) + 1
end
pr("occurrence multiplicity histogram (m occurrences -> number of string pairs): ",
   join(["m=$(p[1])x$(p[2])" for p in sort(collect(multiplicities))], ",  "))
for (key, lst) in sort(collect(occ); by=p -> -length(p[2]))
    pr(@sprintf("  string pair L=%s R=%s : %d occurrences, sum(c_I c_J) = %+.8f  dets %s",
                bits_string(key[1], n_act), bits_string(key[2], n_act), length(lst),
                sum(x[3] for x in lst),
                join(["($(x[1]),$(x[2]))" for x in lst], " ")))
end

# ------------------------------------------------------ 4. assembled matrices
pr("")
pr("# ---------------------------------------------------------------")
pr("# 4. assembled raw NOCI matrices (AHPC production path)")
pr("# ---------------------------------------------------------------")
telemetry = NOCI.ExecutionTelemetry("ahpc-noci", "demo")
mat = NOCI.assemble_ahpc_spin_string_matrices(bfs, one_el, factors, screening;
    telemetry=telemetry)
labels = [b.label for b in bfs]
pr("MEBF order: ", join(labels, ", "))
prmat("S (overlap):", mat.overlap)
prmat("H1e+Enn:", mat.one_body_with_nuclear_repulsion)
prmat("H2e:", mat.two_body)
prmat("H = H1e+Enn+H2e:", mat.hamiltonian)
pr("telemetry counters:")
for (k, v) in sort(collect(telemetry.counters))
    pr("  ", rpad(k, 46), v)
end

# --------------------------------------------------------- 5. solution/couplings
pr("")
pr("# ---------------------------------------------------------------")
pr("# 5. generalized-eigenvalue solution and Loewdin couplings")
pr("# ---------------------------------------------------------------")
result = NOCI.evaluate(input)
prmat("normalized S:", result.normalized.overlap)
prmat("normalized H:", result.normalized.hamiltonian)
pr("energies (Ha):")
for (i, e) in enumerate(result.solution.energies)
    pr(@sprintf("   E%d = %.12f   relative = %.6f eV", i, e, 27.211386245988 * (e - result.solution.energies[1])))
end
prmat("NOCI eigenvectors (columns = states):", result.solution.vectors; digits=6, width=13)
Vp = dataset_pairwise_couplings(result.normalized.hamiltonian, result.normalized.overlap)
Vf = dataset_full_lowdin_couplings(result.normalized.hamiltonian, result.normalized.overlap)
prmat("pairwise symmetric-orthogonalized couplings V_munu (Ha):", Vp; digits=8)
prmat("pairwise couplings (meV):", Vp .* 27211.386245988; digits=3, width=12)
prmat("full-manifold Lowdin H~ off-diagonal (meV):", Vf .* 27211.386245988; digits=3, width=12)

# ------------------------------------------------------ 6. reference comparison
pr("")
pr("# ---------------------------------------------------------------")
pr("# 6. cross-check against archived reference snapshot")
pr("# ---------------------------------------------------------------")
ref = JLD2.load(REF_JLD2)
for (name, mine, key) in (("overlap", mat.overlap, "raw_sbase"),
                          ("one-body", mat.one_body_with_nuclear_repulsion, "raw_h1base"),
                          ("two-body", mat.two_body, "raw_h2base"),
                          ("hamiltonian", mat.hamiltonian, "raw_hbase"))
    pr(@sprintf("  max|%s - ref[%s]| = %.3e", name, key, maximum(abs.(mine .- ref[key]))))
end

# ------------------------------------------------------------------- archive
store["mebf_labels"] = labels
store["raw_matrices"] = (overlap=mat.overlap, one_body=mat.one_body_with_nuclear_repulsion,
                         two_body=mat.two_body, hamiltonian=mat.hamiltonian)
store["normalized_matrices"] = (overlap=result.normalized.overlap, hamiltonian=result.normalized.hamiltonian)
store["energies"] = result.solution.energies
store["couplings_pairwise_ha"] = Vp
store["couplings_full_lowdin_ha"] = Vf
store["basis_overlap"] = one_el.overlap
JLD2.save(joinpath(OUTDIR, "ethene_ahpc_dump.jld2"), store)
pr("")
pr("# dump complete -> ", TXT)
pr("# archive       -> ", joinpath(OUTDIR, "ethene_ahpc_dump.jld2"))
close(io)
println("written: ", TXT)
