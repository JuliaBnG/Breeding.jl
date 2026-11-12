"""
    struct ROCSResult

Solution of a Robust Optimal Contribution Selection (ROCS), Optimal Contribution
Selection (OCS), Robust Truncation Selection (RTS), or Truncation Selection (TS)
problem (Fogg et al., 2026, *Bioinformatics* 42(8):btag569).

# Fields
- `w::Vector{Float64}`: Optimal candidate contribution proportions (`n` candidates).
- `objective::Float64`: Realized objective value ``\\mathbf{w}'\\bar{\\boldsymbol{\\mu}} - \\kappa\\sqrt{\\mathbf{w}'\\boldsymbol{\\Omega}\\mathbf{w}} - \\tfrac{\\lambda}{2}\\mathbf{w}'\\boldsymbol{\\Sigma}\\mathbf{w}``.
- `merit::Float64`: Expected genetic merit ``\\mathbf{w}'\\bar{\\boldsymbol{\\mu}}``.
- `std_error::Float64`: Standard error of expected genetic merit ``\\sqrt{\\mathbf{w}'\\boldsymbol{\\Omega}\\mathbf{w}}``.
- `coancestry::Float64`: Group coancestry ``\\mathbf{w}'\\boldsymbol{\\Sigma}\\mathbf{w}``.
- `λ::Float64`: Penalty weight on group coancestry.
- `κ::Float64`: Robustness parameter penalizing uncertainty in expected genetic merit.
- `iterations::Int`: Number of outer solver iterations executed.
- `converged::Bool`: Whether the KKT optimality tolerance was met.

Supports tuple destructuring `(w, objective, merit, coancestry, std_error) = rocs(...)`.
"""
struct ROCSResult
    w::Vector{Float64}
    objective::Float64
    merit::Float64
    std_error::Float64
    coancestry::Float64
    λ::Float64
    κ::Float64
    iterations::Int
    converged::Bool
end

function Base.iterate(r::ROCSResult, state::Int = 1)
    state == 1 && return (r.w, 2)
    state == 2 && return (r.objective, 3)
    state == 3 && return (r.merit, 4)
    state == 4 && return (r.coancestry, 5)
    state == 5 && return (r.std_error, 6)
    return nothing
end

# Euclidean projection onto a single capped simplex:
#   min_{w} 1/2 ||w - y||_2^2  s.t.  lb <= w <= ub,  sum(w) == target
function _project_capped_simplex!(
    w::AbstractVector{Float64},
    y::AbstractVector{Float64},
    lb::AbstractVector{Float64},
    ub::AbstractVector{Float64},
    target::Float64,
)
    sum_l = sum(lb)
    sum_u = sum(ub)
    target >= sum_l - 1e-12 && target <= sum_u + 1e-12 ||
        throw(ArgumentError("infeasible bounds: sum(lb)=$sum_l, sum(ub)=$sum_u, target=$target"))
    if abs(target - sum_l) <= 1e-14
        copyto!(w, lb)
        return w
    elseif abs(target - sum_u) <= 1e-14
        copyto!(w, ub)
        return w
    end

    # Find shift τ such that sum(clamp.(y .- τ, lb, ub)) == target
    lo = minimum(y .- ub)
    hi = maximum(y .- lb)
    for _ = 1:64
        mid = 0.5 * (lo + hi)
        s = 0.0
        @inbounds @simd for i in eachindex(y)
            s += clamp(y[i] - mid, lb[i], ub[i])
        end
        if s > target
            lo = mid
        else
            hi = mid
        end
        if hi - lo <= 1e-15 * max(1.0, abs(lo), abs(hi))
            break
        end
    end
    τ = 0.5 * (lo + hi)
    @inbounds @simd for i in eachindex(y)
        w[i] = clamp(y[i] - τ, lb[i], ub[i])
    end
    return w
end

function _project_groups!(
    w::Vector{Float64},
    y::AbstractVector{Float64},
    groups::Vector{Vector{Int}},
    targets::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64},
)
    for (g, tgt) in zip(groups, targets)
        _project_capped_simplex!(view(w, g), view(y, g), view(lb, g), view(ub, g), tgt)
    end
    return w
end

# Exact Truncation Selection (λ == 0, κ == 0) within each group:
#   max_{w} μ'w  s.t.  lb <= w <= ub,  sum_{i in g} w_i == target_g
function _solve_ts(
    μ::Vector{Float64},
    groups::Vector{Vector{Int}},
    targets::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64},
)
    n = length(μ)
    w = copy(lb)
    for (g, tgt) in zip(groups, targets)
        rem = tgt - sum(view(lb, g))
        rem >= -1e-12 || throw(ArgumentError("sum(lb) exceeds group target $tgt"))
        sum(view(ub, g)) >= tgt - 1e-12 ||
            throw(ArgumentError("sum(ub) is less than group target $tgt"))
        rem = max(rem, 0.0)
        ord = sortperm(view(μ, g); rev = true)
        for idx in ord
            rem <= 1e-15 && break
            i = g[idx]
            cap = ub[i] - lb[i]
            add = min(rem, cap)
            w[i] += add
            rem -= add
        end
    end
    return w
end

# Active-set KKT polish for strictly convex quadratic subproblem:
#   min_{w} 1/2 w'Hw - c'w  s.t.  sum_{i in g} w_i == target_g,  lb <= w <= ub
function _solve_qp_groups!(
    w::Vector{Float64},
    H::Matrix{Float64},
    c::Vector{Float64},
    groups::Vector{Vector{Int}},
    targets::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64};
    tol::Float64 = 1e-11,
    maxiter::Int = 500,
)
    n = length(w)
    ng = length(groups)
    # State of each variable: -1 = fixed at lb, 0 = free, +1 = fixed at ub
    state = zeros(Int8, n)
    for i = 1:n
        if ub[i] - lb[i] <= 1e-14 || w[i] <= lb[i] + 1e-12
            state[i] = -1
            w[i] = lb[i]
        elseif w[i] >= ub[i] - 1e-12
            state[i] = 1
            w[i] = ub[i]
        end
    end

    # Ensure each group with remaining free mass has at least one free variable
    for (g, tgt) in zip(groups, targets)
        fixed_sum = sum(w[i] for i in g if state[i] != 0; init = 0.0)
        if !any(i -> state[i] == 0, g) && abs(tgt - fixed_sum) > 1e-12
            for i in g
                if ub[i] > lb[i] + 1e-14
                    state[i] = 0
                end
            end
        end
    end

    p = zeros(Float64, n)
    g_grad = zeros(Float64, n)

    for _ = 1:maxiter
        F = findall(==(0), state)
        # Check which groups actually have free variables
        active_groups = Int[]
        m_rhs = Float64[]
        for (gi, g) in enumerate(groups)
            rem = targets[gi]
            has_free = false
            for i in g
                if state[i] == 0
                    has_free = true
                else
                    rem -= w[i]
                end
            end
            if has_free
                push!(active_groups, gi)
                push!(m_rhs, rem)
            end
        end

        w_star = copy(w)
        ν = zeros(Float64, ng)

        if !isempty(F)
            nf = length(F)
            nag = length(active_groups)
            H_FF = Symmetric(H[F, F])
            chol = cholesky(H_FF; check = false)
            if !issuccess(chol)
                H_reg = H[F, F] + (1e-10 * max(1.0, tr(H_FF) / nf)) * I
                chol = cholesky(Symmetric(H_reg))
            end

            rhs_F = copy(c[F])
            Fixed = findall(!=(0), state)
            if !isempty(Fixed)
                mul!(rhs_F, view(H, F, Fixed), view(w, Fixed), -1.0, 1.0)
            end

            # Build M_F (nf × nag)
            M_F = zeros(Float64, nf, nag)
            for (col, gi) in enumerate(active_groups)
                grp_set = Set(groups[gi])
                for (row, fi) in enumerate(F)
                    if fi in grp_set
                        M_F[row, col] = 1.0
                    end
                end
            end

            Hinv_rhs = chol \ rhs_F
            Hinv_M = chol \ M_F
            S = Symmetric(M_F' * Hinv_M)
            ν_sub = S \ (M_F' * Hinv_rhs - m_rhs)
            w_F = Hinv_rhs - Hinv_M * ν_sub
            for (k, fi) in enumerate(F)
                w_star[fi] = w_F[k]
            end
            for (col, gi) in enumerate(active_groups)
                ν[gi] = ν_sub[col]
            end
        end

        p .= w_star .- w
        if norm(p, Inf) <= tol
            # Subspace minimizer reached; check bound multipliers
            w .= clamp.(w_star, lb, ub)
            mul!(g_grad, H, w)
            g_grad .-= c
            for (gi, g) in enumerate(groups)
                if !(gi in active_groups)
                    # Multiplier for a group with no free variables
                    # Pick from average negative gradient on the group
                    ν[gi] = -mean(view(g_grad, g))
                end
                for i in g
                    g_grad[i] += ν[gi]
                end
            end

            worst_val = -tol
            worst_idx = 0
            for i = 1:n
                ub[i] - lb[i] <= 1e-14 && continue
                if state[i] == -1 && g_grad[i] < worst_val
                    worst_val = g_grad[i]
                    worst_idx = i
                elseif state[i] == 1 && -g_grad[i] < worst_val
                    worst_val = -g_grad[i]
                    worst_idx = i
                end
            end

            if worst_idx == 0
                return w, true
            end
            state[worst_idx] = 0
        else
            # Step towards w_star up to first blocking bound
            α = 1.0
            block_idx = 0
            block_state = Int8(0)
            for fi in F
                if p[fi] < -1e-15
                    a = (lb[fi] - w[fi]) / p[fi]
                    if a < α
                        α = a
                        block_idx = fi
                        block_state = Int8(-1)
                    end
                elseif p[fi] > 1e-15
                    a = (ub[fi] - w[fi]) / p[fi]
                    if a < α
                        α = a
                        block_idx = fi
                        block_state = Int8(1)
                    end
                end
            end
            w .+= α .* p
            if block_idx != 0 && α < 1.0 - 1e-14
                state[block_idx] = block_state
                w[block_idx] = block_state == -1 ? lb[block_idx] : ub[block_idx]
            else
                # Snap any boundary hits
                for fi in F
                    if w[fi] <= lb[fi] + 1e-13
                        state[fi] = -1
                        w[fi] = lb[fi]
                    elseif w[fi] >= ub[fi] - 1e-13
                        state[fi] = 1
                        w[fi] = ub[fi]
                    end
                end
            end
        end
    end
    return w, false
end

"""
    rocs(
        μ::AbstractVector{<:Real},
        Σ::AbstractMatrix{<:Real};
        Ω::Union{Nothing, AbstractVecOrMat{<:Real}} = nothing,
        sex::Union{Nothing, AbstractVector} = nothing,
        λ::Real = 1.0,
        κ::Real = 0.0,
        lb::Union{Real, AbstractVector{<:Real}} = 0.0,
        ub::Union{Real, AbstractVector{<:Real}} = 1.0,
        tol::Real = 1e-9,
        maxiter::Integer = 200,
    ) -> ROCSResult

Solve the Robust Optimal Contribution Selection (ROCS) problem of Fogg,
Ortiz-Cuadros, Pocrnić, Hall, and Gorjanc (2026, *Bioinformatics* 42(8):btag569):

```math
\\max_{\\mathbf{l} \\le \\mathbf{w} \\le \\mathbf{u}}
\\left(
    \\mathbf{w}'\\bar{\\boldsymbol{\\mu}}
    - \\kappa \\sqrt{\\mathbf{w}'\\boldsymbol{\\Omega}\\mathbf{w}}
    - \\frac{\\lambda}{2} \\mathbf{w}'\\boldsymbol{\\Sigma}\\mathbf{w}
\\right)
\\quad \\text{subject to} \\quad \\mathbf{M}\\mathbf{w} = \\mathbf{m}
```

where:
- ``\\bar{\\boldsymbol{\\mu}} = \\mathbb{E}(\\boldsymbol{\\mu} \\mid \\mathbf{y})`` is the vector of estimated breeding values (EBVs);
- ``\\boldsymbol{\\Omega} = \\operatorname{Var}(\\boldsymbol{\\mu} \\mid \\mathbf{y})`` is the prediction error covariance (PEC) matrix (or a vector of prediction error variances, ``\\operatorname{diag}(\\boldsymbol{\\Omega})``, e.g. from [`pev`](@ref));
- ``\\boldsymbol{\\Sigma}`` is the genetic relationship (coancestry) matrix among selection candidates (e.g. ``\\mathbf{A}`` or ``\\mathbf{G}``).

# Four Selection Regimes (Table 2 of Fogg et al., 2026)
- `λ = 0, κ = 0`: Classical Truncation Selection (**TS**)
- `λ = 0, κ > 0`: Robust Truncation Selection (**RTS**)
- `λ > 0, κ = 0`: Classical Optimal Contribution Selection (**OCS**)
- `λ > 0, κ > 0`: Robust Optimal Contribution Selection (**ROCS**)

# Arguments
- `μ::AbstractVector{<:Real}`: Candidate EBVs or index values of length `n`.
- `Σ::AbstractMatrix{<:Real}`: Symmetric `n × n` genetic relationship matrix.

# Keywords
- `Ω = nothing`: Prediction error covariance matrix (`n × n`) or prediction error variance vector (`n`). Required when `κ > 0`.
- `sex = nothing`: Optional vector of length `n` partitioning candidates into two sexes (e.g. `0`/`1`, `false`/`true`, or `"M"`/`"F"`), enforcing ``\\sum_{i \\in \\mathcal{S}} w_i = \\tfrac{1}{2}`` and ``\\sum_{i \\in \\mathcal{D}} w_i = \\tfrac{1}{2}``. If `nothing`, enforces the single-cohort constraint ``\\sum_{i=1}^n w_i = 1``.
- `λ::Real = 1.0`: Coancestry penalty weight (``\\lambda \\ge 0``).
- `κ::Real = 0.0`: Robustness uncertainty penalty weight (``\\kappa \\ge 0``).
- `lb = 0.0`: Scalar or vector of lower bounds on ``\\mathbf{w}`` (must be ``\\ge 0``).
- `ub = 1.0`: Scalar or vector of upper bounds on ``\\mathbf{w}``.
- `tol::Real = 1e-9`: KKT optimality tolerance.
- `maxiter::Integer = 200`: Maximum outer SQP / variational QP iterations.

# Returns
- [`ROCSResult`](@ref): Struct containing `w`, `objective`, `merit` (``\\mathbf{w}'\\bar{\\boldsymbol{\\mu}}``), `std_error` (``\\sqrt{\\mathbf{w}'\\boldsymbol{\\Omega}\\mathbf{w}}``), `coancestry` (``\\mathbf{w}'\\boldsymbol{\\Sigma}\\mathbf{w}``), `λ`, `κ`, `iterations`, and `converged`.
"""
function rocs(
    μ_in::AbstractVector{<:Real},
    Σ_in::AbstractMatrix{<:Real};
    Ω::Union{Nothing,AbstractVecOrMat{<:Real}} = nothing,
    sex::Union{Nothing,AbstractVector} = nothing,
    λ::Real = 1.0,
    κ::Real = 0.0,
    lb::Union{Real,AbstractVector{<:Real}} = 0.0,
    ub::Union{Real,AbstractVector{<:Real}} = 1.0,
    tol::Real = 1e-9,
    maxiter::Integer = 200,
)
    n = length(μ_in)
    n >= 1 || throw(ArgumentError("μ must be non-empty"))
    size(Σ_in) == (n, n) ||
        throw(DimensionMismatch("Σ must have size ($n, $n), got $(size(Σ_in))"))
    λ >= 0 || throw(ArgumentError("λ must be non-negative, got $λ"))
    κ >= 0 || throw(ArgumentError("κ must be non-negative, got $κ"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    maxiter >= 1 || throw(ArgumentError("maxiter must be positive"))

    μ = Float64.(μ_in)
    Σ = Matrix{Float64}((Σ_in .+ Σ_in') .* 0.5)

    Ω_mat = if Ω === nothing
        κ > 0 && throw(ArgumentError("Ω must be provided when κ > 0"))
        nothing
    elseif Ω isa AbstractVector
        length(Ω) == n || throw(DimensionMismatch("length(Ω) must be $n"))
        all(>=(0), Ω) || throw(ArgumentError("diagonal entries of Ω must be non-negative"))
         Matrix{Float64}( Diagonal(Float64.(Ω)) )
    else
        size(Ω) == (n, n) || throw(DimensionMismatch("Ω must have size ($n, $n)"))
        Matrix{Float64}((Ω .+ Ω') .* 0.5)
    end

    lb_vec = lb isa Real ? fill(Float64(lb), n) : Float64.(lb)
    ub_vec = ub isa Real ? fill(Float64(ub), n) : Float64.(ub)
    length(lb_vec) == n && length(ub_vec) == n ||
        throw(DimensionMismatch("lb and ub must have length $n"))
    all(lb_vec .>= 0.0) || throw(ArgumentError("lower bounds lb must be non-negative"))
    all(lb_vec .<= ub_vec) || throw(ArgumentError("lb must not exceed ub"))

    groups, targets = if sex === nothing
        ([collect(1:n)], [1.0])
    else
        length(sex) == n || throw(DimensionMismatch("sex must have length $n"))
        u_sex = unique(sex)
        length(u_sex) == 2 ||
            throw(ArgumentError("sex must contain exactly 2 distinct classes, got $(length(u_sex))"))
        g1 = findall(==(u_sex[1]), sex)
        g2 = findall(==(u_sex[2]), sex)
        ([g1, g2], [0.5, 0.5])
    end

    λf = Float64(λ)
    κf = Float64(κ)

    # Helper to evaluate metrics and objective
    eval_metrics(w_vec) = begin
        merit = dot(w_vec, μ)
        coanc = dot(w_vec, Σ * w_vec)
        se = Ω_mat === nothing ? 0.0 : sqrt(max(0.0, dot(w_vec, Ω_mat * w_vec)))
        obj = merit - κf * se - 0.5 * λf * coanc
        (obj, merit, se, coanc)
    end

    # 1. Classical Truncation Selection (λ == 0 && κ == 0)
    if λf == 0.0 && κf == 0.0
        w = _solve_ts(μ, groups, targets, lb_vec, ub_vec)
        obj, merit, se, coanc = eval_metrics(w)
        return ROCSResult(w, obj, merit, se, coanc, λf, κf, 1, true)
    end

    # 2. Initialize feasible w (equal contributions within each group, clamped to [lb, ub])
    w = zeros(Float64, n)
    y0 = zeros(Float64, n)
    for (g, tgt) in zip(groups, targets)
        y0[g] .= tgt / length(g)
    end
    _project_groups!(w, y0, groups, targets, lb_vec, ub_vec)

    # Small ridge for PSD matrices
    scale_Σ = tr(Σ) / n
    scale_Ω = Ω_mat === nothing ? 0.0 : tr(Ω_mat) / n
    ridge = 1e-10 * max(1.0, λf * scale_Σ + κf * scale_Ω)

    # 3. Classical OCS (λ > 0 && κ == 0): single strictly convex QP
    if κf == 0.0
        H = λf .* Σ
        for i = 1:n
            H[i, i] += ridge
        end
        # Warm start with a few spectral projected gradient steps
        grad = similar(w)
        w_trial = similar(w)
        for _ = 1:15
            mul!(grad, H, w)
            grad .-= μ # gradient of 1/2 w'Hw - μ'w
            step = 1.0 / max(1e-8, opnorm(H, 1))
            _project_groups!(w_trial, w .- step .* grad, groups, targets, lb_vec, ub_vec)
            w .= w_trial
        end
        w, conv = _solve_qp_groups!(w, H, μ, groups, targets, lb_vec, ub_vec; tol = Float64(tol))
        obj, merit, se, coanc = eval_metrics(w)
        return ROCSResult(w, obj, merit, se, coanc, λf, κf, 1, conv)
    end

    # 4. Robust Selection (κ > 0: RTS when λ == 0, ROCS when λ > 0)
    # Iteratively solves the exact variational quadratic majorizer / SQP subproblem:
    #   H_k = λ Σ + (κ / sqrt(w_k' Ω w_k)) Ω
    # whose KKT conditions at fixed point w_{k+1} = w_k are identical to Eq. (8) of Fogg et al. (2026).
    H = Matrix{Float64}(undef, n, n)
    w_prev = similar(w)
    d_step = similar(w)
    grad = similar(w)
    w_proj = similar(w)
    conv = false
    iters = 0

    for it = 1:Int(maxiter)
        iters = it
        copyto!(w_prev, w)
        q = max(1e-24, dot(w, Ω_mat * w))
        η = sqrt(q)

        # Assemble effective Hessian H_k = λ Σ + (κ / η) Ω + ridge * I
        α_Ω = κf / η
        @inbounds for j = 1:n, i = 1:n
            H[i, j] = λf * Σ[i, j] + α_Ω * Ω_mat[i, j]
        end
        @inbounds for i = 1:n
            H[i, i] += ridge
        end

        w, qp_conv = _solve_qp_groups!(
            w,
            H,
            μ,
            groups,
            targets,
            lb_vec,
            ub_vec;
            tol = max(1e-12, 0.1 * Float64(tol)),
        )

        # Check KKT stationarity on true objective F(w)
        q_new = max(1e-24, dot(w, Ω_mat * w))
        η_new = sqrt(q_new)
        grad .= μ .- (λf .* (Σ * w)) .- ((κf / η_new) .* (Ω_mat * w))
        _project_groups!(w_proj, w .+ grad, groups, targets, lb_vec, ub_vec)
        kkt_res = norm(w_proj .- w, Inf)
        step_res = norm(w .- w_prev, Inf)

        if kkt_res <= Float64(tol) || (qp_conv && step_res <= Float64(tol))
            conv = true
            break
        end

        # Exact 1D concave line search along d = w - w_prev if needed
        d_step .= w .- w_prev
        f_prev = first(eval_metrics(w_prev))
        f_curr = first(eval_metrics(w))
        if f_curr < f_prev
            α = 0.5
            while α > 1e-6
                w .= w_prev .+ α .* d_step
                first(eval_metrics(w)) >= f_prev && break
                α *= 0.5
            end
        end
    end

    obj, merit, se, coanc = eval_metrics(w)
    return ROCSResult(w, obj, merit, se, coanc, λf, κf, iters, conv)
end

const robust_ocs = rocs
