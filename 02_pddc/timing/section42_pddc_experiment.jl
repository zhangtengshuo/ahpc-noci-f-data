using Dates
using LinearAlgebra
using NOCI
using Printf
using SHA
using Statistics

const ACCURACY_Q_VALUES = (8, 12, 16, 20)
const TIMING_Q_VALUES = (2, 4, 8, 12, 16, 20, 24)
const CONTROLLED_Q = 6
const BIGFLOAT_PRECISION = 256

"""
    deterministic_response_kernel(q)

Construct a deterministic symmetric ordered-response kernel for an occupied rank `q`.
The synthetic values are dimensionless algebraic test data and do not represent molecular electron-repulsion integrals.
"""
function deterministic_response_kernel(q::Int)
    n_response = q * q
    response = Matrix{Float64}(undef, n_response, n_response)
    @inbounds for second in 1:n_response, first in 1:second
        numerator = mod(17 * first + 29 * second + 7 * first * second, 101) - 50
        value = Float64(numerator) / Float64(125 * q)
        response[first, second] = value
        response[second, first] = value
    end
    return response
end

"""
    controlled_factors(q)

Return deterministic, exactly dyadic, unit-triangular factors for the controlled singularity scan.
The factors are well conditioned and all represented entries are exact in binary floating-point arithmetic.
"""
function controlled_factors(q::Int)
    lower = Matrix{Float64}(I, q, q)
    upper = Matrix{Float64}(I, q, q)
    @inbounds for row in 2:q, column in 1:(row - 1)
        sign = isodd(row + 2 * column) ? -1.0 : 1.0
        lower[row, column] = ldexp(sign, -(3 + mod(row + column, 3)))
    end
    @inbounds for column in 2:q, row in 1:(column - 1)
        sign = isodd(2 * row + column) ? -1.0 : 1.0
        upper[row, column] = ldexp(sign, -(3 + mod(row + column + 1, 3)))
    end
    return lower, upper
end

"""
    controlled_samples()

Build the eight q=6 matrices used for the high-precision scan, spanning regular, near-singular, exact-nullity-one, and exact-nullity-two cases.
Each tuple records the sample name, intended nullity, smallest imposed diagonal scale, and matrix.
"""
function controlled_samples()
    q = CONTROLLED_Q
    lower, upper = controlled_factors(q)
    base = Float64[1.0, 1.25, 0.75, 1.125, 0.875, 1.375]
    specifications = [
        ("well_conditioned", 0, 1.0, base[end - 1], base[end]),
        ("near_nullity1_2m8", 0, ldexp(1.0, -8), base[end - 1], ldexp(1.0, -8)),
        ("near_nullity1_2m24", 0, ldexp(1.0, -24), base[end - 1], ldexp(1.0, -24)),
        ("near_nullity1_2m48", 0, ldexp(1.0, -48), base[end - 1], ldexp(1.0, -48)),
        ("exact_nullity1", 1, 0.0, base[end - 1], 0.0),
        ("near_nullity2_2m8_2m12", 0, ldexp(1.0, -12), ldexp(1.0, -8), ldexp(1.0, -12)),
        ("near_nullity2_2m24_2m36", 0, ldexp(1.0, -36), ldexp(1.0, -24), ldexp(1.0, -36)),
        ("exact_nullity2", 2, 0.0, 0.0, 0.0),
    ]
    samples = NamedTuple[]
    for (name, nullity, epsilon, penultimate, final) in specifications
        diagonal = copy(base)
        diagonal[end - 1] = penultimate
        diagonal[end] = final
        matrix = lower * Diagonal(diagonal) * upper
        push!(samples, (; name, nullity, epsilon, matrix))
    end
    return samples
end

"""
    permutation_sign(permutation)

Return the parity sign of a one-based permutation by explicit inversion counting.
This helper is used only by the factorial-cost independent reference.
"""
function permutation_sign(permutation::Vector{Int})
    inversions = 0
    q = length(permutation)
    @inbounds for left in 1:(q - 1), right in (left + 1):q
        inversions += permutation[left] > permutation[right]
    end
    return isodd(inversions) ? -1 : 1
end

"""
    accumulate_explicit_reference!(accumulator, matrix, response, permutation)

Accumulate one permutation's determinant, first derivative, and response-contracted second derivative into a BigFloat reference.
The implementation forms every required product explicitly and never calls a determinant, inverse, production polynomial recurrence, or Hessian-pair routine.
"""
function accumulate_explicit_reference!(accumulator, matrix, response, permutation)
    q = length(permutation)
    sign = BigFloat(permutation_sign(permutation))
    determinant_term = sign
    @inbounds for row in 1:q
        determinant_term *= matrix[row, permutation[row]]
    end
    accumulator.determinant[] += determinant_term

    @inbounds for first_row in 1:q
        first_column = permutation[first_row]
        first_term = sign
        for row in 1:q
            row == first_row && continue
            first_term *= matrix[row, permutation[row]]
        end
        accumulator.first[first_row, first_column] += first_term
        first_response = first_row + (first_column - 1) * q
        for second_row in 1:q
            second_row == first_row && continue
            second_column = permutation[second_row]
            second_response = second_row + (second_column - 1) * q
            second_term = sign
            for row in 1:q
                (row == first_row || row == second_row) && continue
                second_term *= matrix[row, permutation[row]]
            end
            accumulator.second[] += response[first_response, second_response] * second_term
        end
    end
    return nothing
end

"""
    enumerate_reference_permutations!(accumulator, matrix, response, permutation, used, row)

Enumerate all column permutations recursively for the independent BigFloat reference.
The recursion owns only a q-element permutation and Boolean occupancy vector and is restricted to the q=6 validation cases.
"""
function enumerate_reference_permutations!(accumulator, matrix, response, permutation, used, row)
    q = length(permutation)
    if row > q
        accumulate_explicit_reference!(accumulator, matrix, response, permutation)
        return nothing
    end
    @inbounds for column in 1:q
        used[column] && continue
        permutation[row] = column
        used[column] = true
        enumerate_reference_permutations!(accumulator, matrix, response, permutation, used, row + 1)
        used[column] = false
    end
    return nothing
end

"""
    explicit_bigfloat_reference(matrix, response)

Compute q=6 determinant, complete first cofactors, and the contracted second determinant derivative by explicit permutation enumeration at 256-bit precision.
This factorial-cost oracle is independent of every NOCI complementary-cofactor implementation.
"""
function explicit_bigfloat_reference(matrix::Matrix{Float64}, response::Matrix{Float64})
    q = size(matrix, 1)
    q == CONTROLLED_Q || throw(ArgumentError("the explicit reference is restricted to q=$CONTROLLED_Q"))
    return setprecision(BIGFLOAT_PRECISION) do
        matrix_big = BigFloat.(matrix)
        response_big = BigFloat.(response)
        accumulator = (
            determinant=Ref(BigFloat(0)),
            first=zeros(BigFloat, q, q),
            second=Ref(BigFloat(0)),
        )
        enumerate_reference_permutations!(
            accumulator,
            matrix_big,
            response_big,
            zeros(Int, q),
            falses(q),
            1,
        )
        return (
            cofactor_scalar=accumulator.determinant[],
            first_cofactor=accumulator.first,
            contracted_second=accumulator.second[],
        )
    end
end

"""
    diagonal_bigfloat_reference(diagonal, response)

Evaluate the determinant, diagonal first cofactors, and analytic response-contracted second derivative of a diagonal matrix at 256-bit precision.
The formula uses direct products over retained diagonal entries and does not invoke NOCI code.
"""
function diagonal_bigfloat_reference(diagonal::Vector{Float64}, response::Matrix{Float64})
    q = length(diagonal)
    return setprecision(BIGFLOAT_PRECISION) do
        d = BigFloat.(diagonal)
        r = BigFloat.(response)
        determinant = prod(d)
        first = zeros(BigFloat, q, q)
        for index in 1:q
            first[index, index] = prod(d[k] for k in 1:q if k != index)
        end
        second = BigFloat(0)
        for i in 1:q, j in 1:q
            i == j && continue
            response_ii = i + (i - 1) * q
            response_jj = j + (j - 1) * q
            response_ij = i + (j - 1) * q
            response_ji = j + (i - 1) * q
            minor_product = prod(d[k] for k in 1:q if k != i && k != j)
            second += (r[response_ii, response_jj] - r[response_ij, response_ji]) * minor_product
        end
        return (cofactor_scalar=determinant, first_cofactor=first, contracted_second=second)
    end
end

"""
    hessian_pair_context(q, response)

Construct the production state-pair Hessian response and identity union pair maps corresponding to a full q by q ordered response grid.
Response construction is intentionally outside all timed regions.
"""
function hessian_pair_context(q::Int, response::Matrix{Float64})
    n_response = q * q
    coefficients = NOCI.ActiveHamiltonianCoefficients(
        0.0,
        zeros(Float64, n_response),
        0.0,
        zeros(Float64, n_response),
        response,
        zeros(Float64, n_response, n_response),
    )
    state_pair_response = NOCI.state_pair_hessian_response(coefficients, q)
    response_indices = collect(1:n_response)
    pair_count = NOCI._hessian_pair_count(q)
    row_pair_indices = Vector{Int}(undef, pair_count)
    column_pair_indices = Vector{Int}(undef, pair_count)
    NOCI._fill_union_hessian_pair_indices!(row_pair_indices, column_pair_indices, response_indices, q)
    return (; state_pair_response, row_pair_indices, column_pair_indices)
end

"""
    production_hessian_result!(workspace, matrix, context)

Evaluate the production CPU PDDC path from one selected overlap matrix and a preconstructed state-pair Hessian response.
The returned first-cofactor matrix aliases the reusable workspace.
"""
function production_hessian_result!(workspace, matrix, context)
    return NOCI._polynomial_state_pair_hessian_contracted_complementary_cofactors!(
        workspace,
        matrix,
        context.state_pair_response,
        context.row_pair_indices,
        context.column_pair_indices,
    )
end

"""
    error_metrics(candidate, reference)

Return absolute determinant, maximum and Frobenius first-cofactor, and signed/absolute contracted-second errors against a BigFloat reference.
All subtraction is performed in BigFloat before scalar error values are converted to Float64 for CSV output.
"""
function error_metrics(candidate, reference)
    determinant_error = Float64(abs(BigFloat(candidate.cofactor_scalar) - reference.cofactor_scalar))
    first_difference = BigFloat.(candidate.first_cofactor) .- reference.first_cofactor
    first_max_error = Float64(maximum(abs, first_difference))
    first_fro_error = Float64(sqrt(sum(abs2, first_difference)))
    second_signed_error = Float64(BigFloat(candidate.contracted_second) - reference.contracted_second)
    second_abs_error = abs(second_signed_error)
    return (; determinant_error, first_max_error, first_fro_error, second_signed_error, second_abs_error)
end

"""
    validation_pass(candidate, reference, errors; second_rtol=2e-8)

Apply scale-aware numerical acceptance limits to one Float64 PDDC result.
The limits are intentionally looser than observed rounding error and cover exact-rank-deficient references without dividing by a zero target.
"""
function validation_pass(candidate, reference, errors; second_rtol::Float64=2.0e-8)
    determinant_scale = max(1.0, Float64(abs(reference.cofactor_scalar)))
    first_scale = max(1.0, Float64(sqrt(sum(abs2, reference.first_cofactor))))
    second_scale = max(1.0, Float64(abs(reference.contracted_second)))
    return errors.determinant_error <= 2.0e-9 * determinant_scale &&
           errors.first_fro_error <= 1.0e-8 * first_scale &&
           errors.second_abs_error <= second_rtol * second_scale &&
           all(isfinite, candidate.first_cofactor) &&
           isfinite(candidate.cofactor_scalar) &&
           isfinite(candidate.contracted_second)
end

"""
    accuracy_header()

Return the ordered column names for `accuracy.csv`.
"""
function accuracy_header()
    return [
        "section", "sample", "q", "intended_nullity", "epsilon", "sigma_min_over_sigma_max", "reference_kind",
        "c0_reference", "c0_candidate", "c0_abs_error", "c1_max_abs_error", "c1_fro_error",
        "r2_reference", "r2_generic", "r2_hessian", "r2_generic_signed_error", "r2_generic_abs_error",
        "r2_hessian_signed_error", "r2_hessian_abs_error", "generic_hessian_abs_difference", "passed",
    ]
end

"""
    csv_value(value)

Format one machine-readable CSV field without locale-dependent numeric output.
The experiment's labels contain no commas or quotes, so string fields require no additional escaping.
"""
function csv_value(value)
    value === missing && return ""
    value isa Bool && return value ? "true" : "false"
    value isa Integer && return string(value)
    value isa AbstractFloat && return @sprintf("%.17e", Float64(value))
    value isa BigFloat && return @sprintf("%.17e", Float64(value))
    return string(value)
end

"""
    write_csv(path, header, rows)

Write a rectangular collection of named tuples to a CSV file using the supplied ordered header.
Parent directories must already exist and existing files are replaced atomically at the file level.
"""
function write_csv(path::String, header::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, Symbol(column))) for column in header), ','))
        end
    end
    return path
end

"""
    controlled_accuracy_rows()

Run the q=6 controlled scan against the factorial BigFloat oracle and close the generic ordered-response and production Hessian-pair representations.
The returned rows contain raw values and errors for every sample.
"""
function controlled_accuracy_rows()
    q = CONTROLLED_Q
    response = deterministic_response_kernel(q)
    hessian_context = hessian_pair_context(q, response)
    generic_workspace = NOCI.PolynomialContractedCofactorWorkspace{Float64}()
    hessian_workspace = NOCI.PolynomialContractedCofactorWorkspace{Float64}()
    rows = NamedTuple[]
    for sample in controlled_samples()
        reference = explicit_bigfloat_reference(sample.matrix, response)
        generic = NOCI._polynomial_contracted_complementary_cofactors!(generic_workspace, sample.matrix, response)
        generic_first = copy(generic.first_cofactor)
        generic_snapshot = (
            cofactor_scalar=generic.cofactor_scalar,
            first_cofactor=generic_first,
            contracted_second=generic.contracted_second,
        )
        hessian = production_hessian_result!(hessian_workspace, sample.matrix, hessian_context)
        generic_errors = error_metrics(generic_snapshot, reference)
        hessian_errors = error_metrics(hessian, reference)
        singular_values = svdvals(sample.matrix)
        sigma_ratio = minimum(singular_values) / maximum(singular_values)
        closure_error = abs(generic_snapshot.contracted_second - hessian.contracted_second)
        closure_scale = max(1.0, abs(generic_snapshot.contracted_second))
        passed = validation_pass(generic_snapshot, reference, generic_errors) &&
                 validation_pass(hessian, reference, hessian_errors) &&
                 closure_error <= 2.0e-10 * closure_scale
        push!(rows, (
            section="controlled_q6", sample=sample.name, q, intended_nullity=sample.nullity,
            epsilon=sample.epsilon, sigma_min_over_sigma_max=sigma_ratio, reference_kind="explicit_permutation_bigfloat_256",
            c0_reference=reference.cofactor_scalar, c0_candidate=hessian.cofactor_scalar,
            c0_abs_error=hessian_errors.determinant_error, c1_max_abs_error=hessian_errors.first_max_error,
            c1_fro_error=hessian_errors.first_fro_error, r2_reference=reference.contracted_second,
            r2_generic=generic_snapshot.contracted_second, r2_hessian=hessian.contracted_second,
            r2_generic_signed_error=generic_errors.second_signed_error, r2_generic_abs_error=generic_errors.second_abs_error,
            r2_hessian_signed_error=hessian_errors.second_signed_error, r2_hessian_abs_error=hessian_errors.second_abs_error,
            generic_hessian_abs_difference=closure_error, passed,
        ))
    end
    return rows
end

"""
    diagonal_case_diagonals(q)

Return four deterministic diagonal vectors for the large-rank analytic gate: regular, near-singular, exact-nullity-one, and exact-nullity-two.
"""
function diagonal_case_diagonals(q::Int)
    base = [0.875 + Float64(index) / Float64(8 * q) for index in 1:q]
    near = copy(base)
    near[end] = ldexp(1.0, -32)
    nullity_one = copy(base)
    nullity_one[end] = 0.0
    nullity_two = copy(base)
    nullity_two[end - 1:end] .= 0.0
    return [
        (name="well_conditioned", nullity=0, epsilon=minimum(base), diagonal=base),
        (name="near_nullity1_2m32", nullity=0, epsilon=ldexp(1.0, -32), diagonal=near),
        (name="exact_nullity1", nullity=1, epsilon=0.0, diagonal=nullity_one),
        (name="exact_nullity2", nullity=2, epsilon=0.0, diagonal=nullity_two),
    ]
end

"""
    diagonal_accuracy_rows()

Run the q=8,12,16,20 analytic diagonal gates and the required q=20 production Hessian-pair closure.
General ordered-response PDDC is evaluated at every q, while Hessian-pair output is recorded at q=20.
"""
function diagonal_accuracy_rows()
    rows = NamedTuple[]
    for q in ACCURACY_Q_VALUES
        response = deterministic_response_kernel(q)
        generic_workspace = NOCI.PolynomialContractedCofactorWorkspace{Float64}()
        hessian_context = q == 20 ? hessian_pair_context(q, response) : nothing
        hessian_workspace = q == 20 ? NOCI.PolynomialContractedCofactorWorkspace{Float64}() : nothing
        for sample in diagonal_case_diagonals(q)
            matrix = Matrix(Diagonal(sample.diagonal))
            reference = diagonal_bigfloat_reference(sample.diagonal, response)
            generic = NOCI._polynomial_contracted_complementary_cofactors!(generic_workspace, matrix, response)
            generic_snapshot = (
                cofactor_scalar=generic.cofactor_scalar,
                first_cofactor=copy(generic.first_cofactor),
                contracted_second=generic.contracted_second,
            )
            generic_errors = error_metrics(generic_snapshot, reference)
            hessian = q == 20 ? production_hessian_result!(hessian_workspace, matrix, hessian_context) : nothing
            hessian_errors = q == 20 ? error_metrics(hessian, reference) : nothing
            closure_error = q == 20 ? abs(generic_snapshot.contracted_second - hessian.contracted_second) : missing
            closure_scale = max(1.0, abs(generic_snapshot.contracted_second))
            passed = validation_pass(generic_snapshot, reference, generic_errors; second_rtol=1.0e-7) &&
                     (q != 20 || (
                         validation_pass(hessian, reference, hessian_errors; second_rtol=1.0e-7) &&
                         closure_error <= 2.0e-9 * closure_scale
                     ))
            sigma_ratio = minimum(abs, sample.diagonal) / maximum(abs, sample.diagonal)
            push!(rows, (
                section="diagonal_analytic", sample=sample.name, q, intended_nullity=sample.nullity,
                epsilon=sample.epsilon, sigma_min_over_sigma_max=sigma_ratio, reference_kind="diagonal_product_bigfloat_256",
                c0_reference=reference.cofactor_scalar, c0_candidate=q == 20 ? hessian.cofactor_scalar : generic_snapshot.cofactor_scalar,
                c0_abs_error=q == 20 ? hessian_errors.determinant_error : generic_errors.determinant_error,
                c1_max_abs_error=q == 20 ? hessian_errors.first_max_error : generic_errors.first_max_error,
                c1_fro_error=q == 20 ? hessian_errors.first_fro_error : generic_errors.first_fro_error,
                r2_reference=reference.contracted_second, r2_generic=generic_snapshot.contracted_second,
                r2_hessian=q == 20 ? hessian.contracted_second : missing,
                r2_generic_signed_error=generic_errors.second_signed_error, r2_generic_abs_error=generic_errors.second_abs_error,
                r2_hessian_signed_error=q == 20 ? hessian_errors.second_signed_error : missing,
                r2_hessian_abs_error=q == 20 ? hessian_errors.second_abs_error : missing,
                generic_hessian_abs_difference=closure_error, passed,
            ))
        end
    end
    return rows
end

"""
    run_accuracy(output_directory)

Execute every numerical correctness gate, write `accuracy.csv`, and raise if any predeclared acceptance condition fails.
The output combines the q=6 independent permutation reference, all large diagonal gates, and q=20 Hessian-pair closure.
"""
function run_accuracy(output_directory::String)
    rows = vcat(controlled_accuracy_rows(), diagonal_accuracy_rows())
    write_csv(joinpath(output_directory, "accuracy.csv"), accuracy_header(), rows)
    failed = [row for row in rows if !row.passed]
    isempty(failed) || error("$(length(failed)) accuracy cases failed; inspect accuracy.csv")
    return rows
end

"""
    timing_matrices(q, count)

Precompute a deterministic bundle of well-conditioned dense matrices for timing.
Input construction is completed before warm-up and measurement.
"""
function timing_matrices(q::Int, count::Int)
    matrices = Matrix{Float64}[]
    for sample in 1:count
        matrix = Matrix{Float64}(I, q, q)
        @inbounds for column in 1:q, row in 1:q
            matrix[row, column] += 0.006 * sin(Float64(13 * row + 19 * column + 23 * sample)) / sqrt(Float64(q))
        end
        push!(matrices, matrix)
    end
    return matrices
end

"""
    timed_batch!(workspace, matrices, context, task_count)

Measure `task_count` complete production Hessian-pair PDDC calls with a reused workspace and prebuilt inputs.
The returned checksum consumes all three scientific output channels and prevents erroneous call elimination or input reuse assumptions.
"""
function timed_batch!(workspace, matrices, context, task_count::Int)
    checksum = 0.0
    matrix_count = length(matrices)
    start_time = time_ns()
    @inbounds for task in 1:task_count
        matrix = matrices[mod1(task, matrix_count)]
        result = production_hessian_result!(workspace, matrix, context)
        checksum += result.cofactor_scalar + result.first_cofactor[mod1(task, size(matrix, 1))] + result.contracted_second
    end
    elapsed_seconds = Float64(time_ns() - start_time) * 1.0e-9
    return (; elapsed_seconds, checksum)
end

"""
    workspace_payload_bytes(workspace)

Count the exact retained numeric payload in one resized production polynomial workspace by summing the lengths of its arrays.
Julia object headers and caller-owned response data are excluded to match `polynomial_workspace_memory`.
"""
function workspace_payload_bytes(workspace)
    element_count = length(workspace.scaled) + length(workspace.previous) + length(workspace.next) +
                    length(workspace.elementary) + sum(length, workspace.powers) + length(workspace.newton) +
                    length(workspace.projection) + length(workspace.first_cofactor)
    return element_count * sizeof(eltype(workspace.scaled))
end

"""
    run_timing(output_directory)

Measure five rounds for every requested q on one Julia and BLAS thread and write raw `timing.csv` plus disaggregated `memory.csv`.
Each q uses a pilot to choose a common round task count targeting `PDDC_TIMING_TARGET_SECONDS`, defaulting to one second.
"""
function run_timing(output_directory::String)
    Threads.nthreads() == 1 || error("timing requires exactly one Julia thread")
    BLAS.set_num_threads(1)
    target_seconds = parse(Float64, get(ENV, "PDDC_TIMING_TARGET_SECONDS", "1.0"))
    target_seconds > 0.0 || throw(ArgumentError("PDDC_TIMING_TARGET_SECONDS must be positive"))
    timing_rows = NamedTuple[]
    memory_rows = NamedTuple[]
    for q in TIMING_Q_VALUES
        response = deterministic_response_kernel(q)
        context = hessian_pair_context(q, response)
        matrices = timing_matrices(q, 8)
        workspace = NOCI.PolynomialContractedCofactorWorkspace{Float64}()
        for matrix in matrices
            production_hessian_result!(workspace, matrix, context)
        end
        pilot_tasks = max(1, min(32, 256 ÷ q))
        pilot = timed_batch!(workspace, matrices, context, pilot_tasks)
        seconds_per_task = max(pilot.elapsed_seconds / pilot_tasks, eps(Float64))
        task_count = clamp(ceil(Int, target_seconds / seconds_per_task), length(matrices), 2_000_000)
        for round in 1:5
            GC.gc()
            measured = timed_batch!(workspace, matrices, context, task_count)
            push!(timing_rows, (
                q, round, task_count, total_seconds=measured.elapsed_seconds,
                seconds_per_task=measured.elapsed_seconds / task_count, checksum=measured.checksum,
            ))
        end
        estimate = NOCI.polynomial_workspace_memory(Float64, q)
        payload = workspace_payload_bytes(workspace)
        payload == estimate.host_bytes || error("workspace payload mismatch at q=$q: $payload != $(estimate.host_bytes)")
        push!(memory_rows, (
            q, polynomial_workspace_bytes=estimate.host_bytes, measured_numeric_payload_bytes=payload,
            hessian_pair_response_bytes=sizeof(eltype(context.state_pair_response)) * length(context.state_pair_response),
            hessian_pair_count=size(context.state_pair_response, 1),
        ))
    end
    write_csv(
        joinpath(output_directory, "timing.csv"),
        ["q", "round", "task_count", "total_seconds", "seconds_per_task", "checksum"],
        timing_rows,
    )
    write_csv(
        joinpath(output_directory, "memory.csv"),
        ["q", "polynomial_workspace_bytes", "measured_numeric_payload_bytes", "hessian_pair_response_bytes", "hessian_pair_count"],
        memory_rows,
    )
    return timing_rows, memory_rows
end

"""
    fitted_timing_slope(timing_rows)

Fit the observed log-log timing slope over the predeclared high-q interval q=8,12,16,20,24 using per-q round medians.
The statistic characterizes finite-size measurements and is not used as a proof of formal complexity.
"""
function fitted_timing_slope(timing_rows)
    selected_q = [8, 12, 16, 20, 24]
    medians = [median(row.seconds_per_task for row in timing_rows if row.q == q) for q in selected_q]
    x = log.(Float64.(selected_q))
    y = log.(medians)
    centered_x = x .- mean(x)
    slope = dot(centered_x, y .- mean(y)) / dot(centered_x, centered_x)
    return (; selected_q, medians, slope)
end

"""
    cpu_model_name()

Read the first Linux `/proc/cpuinfo` model-name entry, returning `unknown` when the field is unavailable.
"""
function cpu_model_name()
    isfile("/proc/cpuinfo") || return "unknown"
    for line in eachline("/proc/cpuinfo")
        startswith(line, "model name") || continue
        fields = split(line, ':'; limit=2)
        return length(fields) == 2 ? strip(fields[2]) : strip(line)
    end
    return "unknown"
end

"""
    json_escape(value)

Escape a string for the limited flat JSON metadata emitted by this experiment.
"""
function json_escape(value)
    text = string(value)
    text = replace(text, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "\\r", "\t" => "\\t")
    return text
end

"""
    write_run_json(path, mode, accuracy_rows, timing_rows)

Write source, CPU, Julia, BLAS, thread, invocation, and validation metadata to a flat JSON file.
Git identity is injected by the launcher because the deployed GV800 source snapshot intentionally excludes `.git`.
"""
function write_run_json(path::String, mode::String, accuracy_rows, timing_rows)
    slope = isempty(timing_rows) ? missing : fitted_timing_slope(timing_rows).slope
    fields = [
        "timestamp_utc" => string(now(UTC)),
        "mode" => mode,
        "noci_source_commit" => get(ENV, "NOCI_SOURCE_COMMIT", "unknown"),
        "noci_source_status" => get(ENV, "NOCI_SOURCE_STATUS", "unknown"),
        "experiment_source_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "hostname" => gethostname(),
        "cpu_model" => cpu_model_name(),
        "logical_cpu_count" => Sys.CPU_THREADS,
        "julia_version" => string(VERSION),
        "julia_executable" => joinpath(Sys.BINDIR, Base.julia_exename()),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "blas_configuration" => sprint(show, BLAS.get_config()),
        "julia_depot_path" => get(ENV, "JULIA_DEPOT_PATH", ""),
        "julia_pkg_offline" => get(ENV, "JULIA_PKG_OFFLINE", ""),
        "timing_target_seconds" => get(ENV, "PDDC_TIMING_TARGET_SECONDS", ""),
        "command_arguments" => join(ARGS, " "),
        "accuracy_case_count" => length(accuracy_rows),
        "accuracy_pass_count" => count(row -> row.passed, accuracy_rows),
        "timing_round_count" => length(timing_rows),
        "high_q_timing_slope" => slope,
    ]
    open(path, "w") do io
        println(io, "{")
        for (index, (key, value)) in enumerate(fields)
            suffix = index == length(fields) ? "" : ","
            if value === missing
                println(io, "  \"$(json_escape(key))\": null$suffix")
            elseif value isa Number
                println(io, "  \"$(json_escape(key))\": $(value)$suffix")
            else
                println(io, "  \"$(json_escape(key))\": \"$(json_escape(value))\"$suffix")
            end
        end
        println(io, "}")
    end
    return path
end

"""
    main(arguments)

Run `accuracy`, `timing`, or `all` into a required output directory and emit frozen CSV and JSON artifacts.
No package manager API, network operation, source mutation, or runtime package installation is performed.
"""
function main(arguments)
    length(arguments) == 2 || error("usage: section42_pddc_experiment.jl OUTPUT_DIRECTORY accuracy|timing|all")
    output_directory = abspath(arguments[1])
    mode = arguments[2]
    mode in ("accuracy", "timing", "all") || throw(ArgumentError("mode must be accuracy, timing, or all"))
    mkpath(output_directory)
    BLAS.set_num_threads(1)
    accuracy_rows = mode in ("accuracy", "all") ? run_accuracy(output_directory) : NamedTuple[]
    timing_rows, memory_rows = mode in ("timing", "all") ? run_timing(output_directory) : (NamedTuple[], NamedTuple[])
    write_run_json(joinpath(output_directory, "run.json"), mode, accuracy_rows, timing_rows)
    println("PDDC_ACCURACY_CASES=$(length(accuracy_rows))")
    println("PDDC_ACCURACY_PASSED=$(count(row -> row.passed, accuracy_rows))")
    println("PDDC_TIMING_ROUNDS=$(length(timing_rows))")
    if !isempty(timing_rows)
        fitted = fitted_timing_slope(timing_rows)
        println("PDDC_HIGH_Q_TIMING_SLOPE=$(fitted.slope)")
    end
    println("PDDC_OUTPUT_DIRECTORY=$output_directory")
    return nothing
end

main(ARGS)
