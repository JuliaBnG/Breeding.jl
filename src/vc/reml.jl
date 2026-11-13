"""
    ResidualStructure

Residual (co)variance structure for [`reml`](@ref):
[`Homogeneous`](@ref), [`HetClasses`](@ref) or [`MultiTrait`](@ref).
"""
abstract type ResidualStructure end

"""
    Homogeneous(σ²)

Homogeneous residual variance, ``\\operatorname{var}(e) = I\\sigma^2``.
"""
struct Homogeneous <: ResidualStructure
    σ²::Float64
end

"""
    HetClasses(class, σ²)
    HetClasses(t, ranges, σ²)

Heterogeneous residual variances: record `i` has variance
`σ²[class[i]]`, classes `1:length(σ²)`, records independent
(``\\operatorname{var}(e) = \\bigoplus_c I_{n_c}\\sigma^2_c``). The
second form assigns record `i` to the first of `ranges` that contains
`t[i]`, e.g. classes of days in milk for random-regression models
(WOMBAT `VAR residual 1 HET`).
"""
struct HetClasses <: ResidualStructure
    class::Vector{Int}
    σ²::Vector{Float64}
    function HetClasses(class::AbstractVector{<:Integer}, σ²::AbstractVector{<:Real})
        m = length(σ²)
        all(c -> 1 ≤ c ≤ m, class) || throw(ArgumentError("classes must be in 1:$m"))
        all(c -> any(==(c), class), 1:m) || throw(ArgumentError("every class needs records"))
        new(Vector{Int}(class), Vector{Float64}(σ²))
    end
end

function HetClasses(t::AbstractVector, ranges::AbstractVector, σ²::AbstractVector{<:Real})
    class = map(t) do x
        c = findfirst(r -> x in r, ranges)
        c === nothing && throw(ArgumentError("value $x is in none of the ranges"))
        c
    end
    HetClasses(class, σ²)
end

"""
    MultiTrait(trait, group, R0; class = nothing)

Residual covariances among traits, for multi-trait (and multi-trait
random-regression) models with the observations in long format:
observation `i` is a record of trait `trait[i]` ∈ `1:t` in residual
group `group[i]` (any values, e.g. `(animal, time)` for records taken
together). Residuals in the same group have covariance matrix
`R0[S, S]`, `S` the traits recorded in the group (each at most once);
residuals in different groups are independent, so missing records need
no special treatment. `R0` is the ``t \\times t`` starting matrix, or a
vector of them with `class` (per observation, constant within a group)
for heterogeneous residual covariance matrices. Every pair of traits
must be recorded together in some group of every class.

WOMBAT: `VAR residual t` of a multivariate analysis with `RPTCOV
TSELECT` or `ALIGNED` (groups: animal × time), and `VAR residual 1` of a
multivariate random regression (groups: subject × control variable).
"""
struct MultiTrait <: ResidualStructure
    trait::Vector{Int}
    group::Vector{Int}
    class::Vector{Int}
    R0::Vector{Matrix{Float64}}
end

function MultiTrait(trait::AbstractVector{<:Integer}, group::AbstractVector,
                    R0::AbstractVector{<:AbstractMatrix}; class = nothing)
    n = length(trait)
    length(group) == n || throw(DimensionMismatch("$(length(group)) groups for $n observations"))
    t = size(R0[1], 1)
    all(R -> size(R) == (t, t), R0) || throw(DimensionMismatch("every R0 must be $t × $t"))
    all(R -> isposdef(Symmetric(Matrix{Float64}(R))), R0) ||
        throw(ArgumentError("R0 must be positive definite"))
    all(x -> 1 ≤ x ≤ t, trait) || throw(ArgumentError("traits must be in 1:$t"))
    cls = class === nothing ? ones(Int, n) : Vector{Int}(class)
    length(cls) == n || throw(DimensionMismatch("$(length(cls)) classes for $n observations"))
    all(c -> 1 ≤ c ≤ length(R0), cls) || throw(ArgumentError("classes must be in 1:$(length(R0))"))
    ids = Dict{eltype(group),Int}()
    g = [get!(ids, k, length(ids) + 1) for k in group]
    allunique(zip(g, trait)) || throw(ArgumentError("a residual group has a trait more than once"))
    MultiTrait(Vector{Int}(trait), g, cls, [Matrix{Float64}(R) for R in R0])
end

MultiTrait(trait::AbstractVector{<:Integer}, group::AbstractVector, R0::AbstractMatrix; class = nothing) =
    MultiTrait(trait, group, [R0]; class = class)

"""
    MultiTrait(obs::AbstractMatrix{Bool}, R0; class = nothing)

The same for records stacked from an `n × t` table by
[`mt_stack`](@ref) (trait-major): every row of the table is a residual
group, `obs` its observation pattern; `class` (optional) is per row.
"""
function MultiTrait(obs::AbstractMatrix{Bool}, R0; class = nothing)
    trait = [j for j in axes(obs, 2) for i in axes(obs, 1) if obs[i, j]]
    row = [i for j in axes(obs, 2) for i in axes(obs, 1) if obs[i, j]]
    MultiTrait(trait, row, R0 isa AbstractMatrix ? [R0] : R0;
               class = class === nothing ? nothing : class[row])
end

# Residual blocks: the groups of class `c` with recorded traits `P`
# (sorted); `rows[a, g]` is the observation of trait `P[a]` in the g-th
# of these groups. Single-trait structures are the case t = 1.
struct _RBlock
    c::Int
    P::Vector{Int}
    rows::Matrix{Int}
end

_rmats(r::Homogeneous) = [fill(r.σ², 1, 1)]
_rmats(r::HetClasses) = [fill(s, 1, 1) for s in r.σ²]
_rmats(r::MultiTrait) = [copy(R) for R in r.R0]

_rblocks(::Homogeneous, n) = [_RBlock(1, [1], reshape(collect(1:n), 1, n))]

function _rblocks(r::HetClasses, n)
    length(r.class) == n || throw(DimensionMismatch("$(length(r.class)) classes for $n records"))
    [_RBlock(c, [1], reshape(findall(==(c), r.class), 1, :)) for c in eachindex(r.σ²)]
end

function _rblocks(r::MultiTrait, n)
    length(r.trait) == n || throw(DimensionMismatch("$(length(r.trait)) traits for $n records"))
    t = size(r.R0[1], 1)
    members = Dict{Int,Vector{Int}}()
    for i in 1:n
        push!(get!(members, r.group[i], Int[]), i)
    end
    cols = Dict{Tuple{Int,Vector{Int}},Vector{Vector{Int}}}()
    for g in sort!(collect(keys(members)))
        obs = sort!(members[g]; by = i -> r.trait[i])
        P = r.trait[obs]
        allunique(P) || throw(ArgumentError("a residual group has a trait more than once"))
        c = r.class[obs[1]]
        all(==(c), r.class[obs]) || throw(ArgumentError("class not constant within a residual group"))
        push!(get!(cols, (c, P), Vector{Int}[]), obs)
    end
    blocks = [_RBlock(c, P, reshape(reduce(vcat, v), length(P), :)) for ((c, P), v) in sort!(collect(cols); by = first)]
    for c in eachindex(r.R0)
        seen = falses(t, t)
        for b in blocks
            b.c == c && (seen[b.P, b.P] .= true)
        end
        if !all(seen)
            bad = [(i, j) for j in 1:t for i in j:t if !seen[i, j]]
            throw(ArgumentError("class $c: traits never recorded together in a residual group: $bad"))
        end
    end
    blocks
end

"""
    REMLResult

Output of [`reml`](@ref), all at the returned estimates:

- `σ²e`: the residual variance (homogeneous residuals; `NaN`
  otherwise); `R0`: the residual covariance matrix of every class
  (``1 \\times 1`` for single-trait structures); `residual`: the
  residual parameters (the lower triangle of each `R0` by columns, i.e.
  the variance of every class for [`HetClasses`](@ref));
- `G0`: the covariance matrices of the random effects;
- `logL`: the restricted log-likelihood (Mrode's form, without
  constants, as in WOMBAT);
- `AIinv`: the inverse of the average-information matrix on the
  covariance scale (approximate sampling covariances); parameter order:
  `residual`, then the lower triangle of each `G0`, by columns;
- `history`: parameters per iterate (one column each, starting with the
  initial values), `logLs`: their log-likelihoods;
- `iterations`, and the mixed-model `solutions`;
- `steps`: for every iterate, the fraction of the AI step taken (1 for
  a full step, ``2^{-h}`` after ``h`` halvings), or 0 for an EM step.
"""
struct REMLResult
    σ²e::Float64
    G0::Vector{Matrix{Float64}}
    logL::Float64
    AIinv::Matrix{Float64}
    history::Matrix{Float64}
    logLs::Vector{Float64}
    iterations::Int
    solutions::Vector{Float64}
    residual::Vector{Float64}
    steps::Vector{Float64}
    R0::Vector{Matrix{Float64}}
end

# ---- parameters ------------------------------------------------------
# θ: the lower triangles (by columns) of all covariance matrices,
# residual matrices first, then the G₀ₖ; `dims` their orders. On the
# Cholesky scale, 1 × 1 residual matrices are parameterised by log σ²
# (`logv`), all others by their Cholesky factor with log diagonal.

# lower-triangle index pairs of an m × m covariance matrix
_vech(m) = [(i, j) for j in 1:m for i in j:m]

_pars(mats) = reduce(vcat, [[M[i, j] for (i, j) in _vech(size(M, 1))] for M in mats]; init = Float64[])

function _unpack(θ, dims)
    mats = Matrix{Float64}[]
    k = 0
    for m in dims
        M = zeros(m, m)
        for (i, j) in _vech(m)
            k += 1
            M[i, j] = M[j, i] = θ[k]
        end
        push!(mats, M)
    end
    mats
end

_valid(θ, dims) = all(M -> isposdef(Symmetric(M)), _unpack(θ, dims))

function _to_phi(θ, dims, logv)
    φ = Float64[]
    for (M, lv) in zip(_unpack(θ, dims), logv)
        if lv
            push!(φ, log(M[1, 1]))
        else
            L = cholesky(Symmetric(M)).L
            append!(φ, [i == j ? log(L[i, i]) : L[i, j] for (i, j) in _vech(size(M, 1))])
        end
    end
    φ
end

function _from_phi(φ, dims, logv)
    θ = Float64[]
    k = 0
    for (m, lv) in zip(dims, logv)
        if lv
            push!(θ, exp(φ[k+=1]))
        else
            L = zeros(m, m)
            for (i, j) in _vech(m)
                k += 1
                L[i, j] = i == j ? exp(φ[k]) : φ[k]
            end
            G = L * L'
            append!(θ, [G[i, j] for (i, j) in _vech(m)])
        end
    end
    θ
end

# keep every Cholesky pivot at least `rel` times the largest one of its
# matrix (WOMBAT's limit for "small pivots"), so that estimates on the
# boundary of the parameter space stay numerically positive definite
function _floor_phi!(φ, dims, logv; rel = 1e-4)
    k = 0
    for (m, lv) in zip(dims, logv)
        pr = _vech(m)
        if !lv
            d = [k + c for (c, (i, j)) in enumerate(pr) if i == j]
            lo = maximum(φ[d]) + log(rel)
            φ[d] .= max.(φ[d], lo)
        end
        k += length(pr)
    end
    φ
end

# Jacobian ∂θ/∂φ (block diagonal)
function _jacobian(θ, dims, logv)
    blocks = Matrix{Float64}[]
    for (M, lv) in zip(_unpack(θ, dims), logv)
        if lv
            push!(blocks, fill(M[1, 1], 1, 1))
            continue
        end
        m = size(M, 1)
        L = cholesky(Symmetric(M)).L
        pr = _vech(m)
        J = zeros(length(pr), length(pr))
        for (c, (i, j)) in enumerate(pr)
            for (r, (a, b)) in enumerate(pr)     # ∂G_ab/∂L_ij = δ_ai L_bj + L_aj δ_bi
                J[r, c] = (a == i ? L[b, j] : 0.0) + (b == i ? L[a, j] : 0.0)
            end
            i == j && (J[:, c] .*= L[i, i])      # log-transformed diagonal
        end
        push!(blocks, J)
    end
    cat(blocks...; dims = (1, 2))
end

"""
    reml(y, X, effects; σ²e | residual, method = :ai, param = :natural,
         zero = Int[], fix = Int[], equal = Vector{Int}[], tol = 1e-8, ltol = nothing,
         ntol = nothing, maxiter = 200, maxhalve = 10, maxstep = 5.0,
         backend = :auto) -> REMLResult

Restricted maximum likelihood (Patterson and Thompson, 1971) for the
mixed model

```math
y = Xb + \\sum_k Z_k u_k + e, \\quad \\operatorname{var}(u_k) = G_{0k} \\otimes K_k,
\\quad \\operatorname{var}(e) = R,
```

where each [`RandomEffect`](@ref) in `effects` supplies ``Z_k``,
``K_k^{-1}`` and the starting ``G_{0k}`` (so correlated direct–maternal
effects or random-regression covariance matrices are estimated too).
The residual structure is given either as `σ²e` (starting value of a
homogeneous variance) or as `residual`, a [`Homogeneous`](@ref),
[`HetClasses`](@ref) or [`MultiTrait`](@ref) with starting values.
`zero` lists constrained fixed effects (as in [`solve_mme`](@ref)).

Multi-trait models use observations in long format (one element of `y`
per record of one trait): trait-specific fixed effects are columns of
`X`, and a random effect covering several traits has one component per
trait (and regression coefficient), its `Z` having nonzeros only for the
records of that trait; ``G_{0k}`` is then the covariance matrix among
all of them. The residual structure is a [`MultiTrait`](@ref).

`method`:
- `:ai` – average information (Gilmour et al., 1995): first derivatives
  from MME solutions and traces of ``C^{-1}``, and
  ``AI_{ij} = \\tfrac12 f_i'Pf_j`` from working variates
  ``f_i = (\\partial V/\\partial\\theta_i)Py`` (``(\\partial R/\\partial\\theta_i)
  R^{-1}\\hat e`` for residual parameters, ``Z_k(\\partial G_{0k}
  G_{0k}^{-1} \\otimes I)\\hat u_k`` for ``G_{0k}``). The AI matrix is
  made safely positive definite by adding ``\\tau I``, ``\\tau = \\max(0,
  \\delta - \\lambda_{min})``, ``\\delta = 3 \\times 10^{-6}\\lambda_{max}``
  (WOMBAT's default modification). A step must not decrease the
  log-likelihood; otherwise it is halved up to `maxhalve` times, then
  replaced by an EM step.
- `:em` – expectation maximization: for homogeneous residuals
  ``(n - p)\\sigma^2_e = e'y``; otherwise, for every residual class,
  ``R_0 = N^{-1}\\sum_g E(e_g e_g' \\mid y)`` over its ``N`` groups, with
  ``E(e_ge_g') = B_g(\\hat e_g\\hat e_g' + W_gC^{-1}W_g')B_g'`` plus the
  conditional covariance of the residuals of traits not recorded in the
  group (``B_g`` predicts them from those recorded; for single-trait
  classes ``n_c\\sigma^2_c = e_c'e_c + \\operatorname{tr}(C^{-1}W_c'W_c)``);
  and ``G_{0k} = (\\hat U_k'K_k^{-1}\\hat U_k + T_k)/q_k`` with ``T_k[i,j] =
  \\operatorname{tr}(C^{k_ik_j}K_k^{-1})``.

Constraints (with `param = :natural`): `fix` lists parameters (indices
into the parameter vector, ordered as `AIinv`) kept at their starting
values, e.g. zero covariances between blocks of a covariance matrix;
`equal` lists sets of parameters estimated as one value (started at
their mean), e.g. the same variance in several classes. Steps are taken
on the free parameters ``\beta`` with ``\theta = \theta_0 + K\beta``,
gradient ``K'g`` and average information ``K'AK``; `AIinv` is then
``K(K'AK)^{-1}K'``. A fall-back EM step is projected onto the
constraints and taken only if it does not decrease the log-likelihood.

`param` (AI only): `:natural` takes AI steps on the covariances
themselves (steps leaving the parameter space are halved);
`:cholesky` on log residual variances (single-trait structures) and
the Cholesky factors of all other covariance matrices with
log-transformed diagonals, which keeps every iterate positive definite.
The average information is computed on the covariance scale and
transformed with the Jacobian (as in WOMBAT). With `:cholesky`, no
element of a step may exceed `maxstep` (the step is scaled down), and
every Cholesky pivot is kept at least ``10^{-4}`` times the largest
pivot of its matrix, so that estimates on the boundary of the parameter
space remain numerically positive definite.

Iteration stops when the largest change in a covariance parameter is
below `tol` times the largest parameter (at least 1); or, if `ltol` is
given, when an AI step increases the log-likelihood by less than `ltol`
(WOMBAT uses ``5 \\times 10^{-4}``) and, if `ntol` is given, the
increase expected from the next AI step, ``\\tfrac12 g'H^{-1}g`` (half the
Newton decrement, with ``g`` the gradient and ``H`` the modified AI
matrix on the scale of `param`), is below `ntol`; or after `maxiter`
iterates. On flat likelihood surfaces a single small increase does not
mean that the maximum is near; `ntol` guards against stopping early.
Near the boundary of the parameter space the likelihood can be very
flat in some directions, and the parameter criterion alone may need
many iterates.

The log-likelihood is ``-\\tfrac12[y'Py + \\log|V| + \\log|X'V^{-1}X|]``
computed through ``\\log|C| + \\log|R| + \\sum_k \\log|G_k|``.

`backend`: both back-ends assemble ``C`` on one sparsity pattern kept
over all iterates. `:dense` then inverts ``C`` densely (small models);
`:sparse` reuses the symbolic CHOLMOD factorisation, takes the traces
from the selected inverse of ``C`` on its own pattern (Takahashi
equations, via SelectedInversion.jl) and the AI matrix from solves with
the factor; `:auto` (default) uses `:sparse` above 2000 equations.
Dense relationship matrices (e.g. ``G^{-1}``) make the corresponding
blocks of ``C`` and its factor dense.
"""
function reml(y::AbstractVector, X::AbstractMatrix, effects::AbstractVector{<:RandomEffect};
              σ²e = nothing, residual = nothing, method::Symbol = :ai, param::Symbol = :natural,
              zero = Int[], fix = Int[], equal = Vector{Int}[], tol = 1e-8, ltol = nothing,
              ntol = nothing, maxiter = 200, maxhalve = 10, maxstep = 5.0, backend::Symbol = :auto)
    method in (:ai, :em) || throw(ArgumentError("method must be :ai or :em"))
    param in (:natural, :cholesky) || throw(ArgumentError("param must be :natural or :cholesky"))
    backend in (:auto, :dense, :sparse) || throw(ArgumentError("backend must be :auto, :dense or :sparse"))
    (σ²e === nothing) == (residual === nothing) &&
        throw(ArgumentError("give exactly one of σ²e or residual"))
    res = residual === nothing ? Homogeneous(σ²e) : residual
    hom = res isa Homogeneous
    n = length(y)
    yv = Vector{Float64}(y)
    blocks = _rblocks(res, n)
    Rm = _rmats(res)
    nR = length(Rm)
    G0 = [copy(r.G0) for r in effects]
    dims = [size.(Rm, 1); size.(G0, 1)]
    logv = [k ≤ nR && dims[k] == 1 for k in eachindex(dims)]
    W = SparseMatrixCSC{Float64,Int}(hcat(sparse(X), (r.Z for r in effects)...))
    # log|K_k| from K_k⁻¹ (constant over iterations)
    logK = [-_logdet_spd(r.Kinv) for r in effects]
    store = _setup(yv, X, W, effects, blocks, zero)
    neq = length(store.keepcol)
    p = neq - sum(size(r.Z, 2) for r in effects; init = 0)
    be = backend == :auto ? (neq > 2000 ? :sparse : :dense) : backend
    θ = _pars([Rm; G0])
    K = _constraints(length(θ), fix, equal)
    if K !== nothing
        param == :natural && method == :ai ||
            throw(ArgumentError("fix and equal need method = :ai and param = :natural"))
        for g in equal
            θ[g] .= sum(θ[g]) / length(g)
        end
    end
    cons = (; K, fix, equal)
    ctx = (; yv, W, effects, blocks, nR, dims, logv, logK, hom, p, be, store, cons)

    ev = _reml_eval(ctx, θ; full = true)
    hist = [θ]
    logLs = [ev.logL]
    steps = Float64[]
    it = 0
    while it < maxiter
        θn, λ = method == :em ? (_pars([ev.em_R; ev.em_G]), 0.0) :
                _ai_step(ctx, θ, ev, param, maxhalve, maxstep)
        evn = _reml_eval(ctx, θn; full = true)
        conv = maximum(abs, θn - θ) < tol * max(1, maximum(abs, θ)) ||
               (ltol !== nothing && λ > 0 && evn.logL - ev.logL < ltol &&
                (ntol === nothing || _gain(evn, θn, ctx, param) < ntol))
        θ, ev = θn, evn
        push!(hist, θ)
        push!(logLs, ev.logL)
        push!(steps, λ)
        it += 1
        conv && break
    end
    mats = _unpack(θ, dims)
    R = mats[1:nR]
    nθR = sum(d * (d + 1) ÷ 2 for d in dims[1:nR])
    AIinv = K === nothing ? inv(ev.AI) : K * (Symmetric(K' * ev.AI * K) \ K')
    REMLResult(hom ? R[1][1, 1] : NaN, mats[nR+1:end], ev.logL, AIinv, reduce(hcat, hist), logLs, it,
               ev.x, θ[1:nθR], steps, R)
end

# θ = θ₀ + Kβ: one column of K per free parameter and per set in `equal`
function _constraints(np, fix, equal)
    isempty(fix) && isempty(equal) && return nothing
    inset = reduce(vcat, [collect(Int, fix); [collect(Int, g) for g in equal]])
    all(i -> 1 ≤ i ≤ np, inset) || throw(ArgumentError("constraint index outside 1:$np"))
    allunique(inset) || throw(ArgumentError("a parameter is in more than one constraint"))
    cols = Vector{Float64}[]
    for i in 1:np
        i in inset && continue
        v = zeros(np); v[i] = 1
        push!(cols, v)
    end
    for g in equal
        v = zeros(np); v[g] .= 1
        push!(cols, v)
    end
    reduce(hcat, cols)
end

# project parameters onto the constraints (fixed at θ₀, equal sets averaged)
function _project(θ, θ0, cons)
    c = copy(θ)
    c[cons.fix] = θ0[cons.fix]
    for g in cons.equal
        c[g] .= sum(c[g]) / length(g)
    end
    c
end

# WOMBAT's default modification of the AI matrix: add τI so that the
# smallest eigenvalue is at least 3e-6 times the largest
function _modify_ai(A)
    λ = eigvals(Symmetric(Matrix(A)))
    τ = max(0.0, 3e-6 * λ[end] - λ[1])
    Symmetric(Matrix(A) + τ * I)
end

# increase in log L expected from a (modified) AI step: ½ g′H⁻¹g
function _gain(ev, θ, ctx, param)
    g, H = ev.grad, ev.AI
    if ctx.cons.K !== nothing
        K = ctx.cons.K
        g, H = K' * g, K' * H * K
    elseif param == :cholesky
        J = _jacobian(θ, ctx.dims, ctx.logv)
        g, H = J' * g, J' * H * J
    end
    dot(g, _modify_ai(H) \ g) / 2
end

# AI step from θ, halved until the log-likelihood does not decrease; EM
# step if no step length does. Returns the new parameters and the
# fraction of the AI step taken (0 for the EM step).
function _ai_step(ctx, θ, ev, param, maxhalve, maxstep)
    (; dims, logv) = ctx
    K = ctx.cons.K
    if K !== nothing
        step = K * (_modify_ai(K' * ev.AI * K) \ (K' * ev.grad))
        cand = λ -> θ + λ * step
        ok = c -> _valid(c, dims)
    elseif param == :natural
        step = _modify_ai(ev.AI) \ ev.grad
        cand = λ -> θ + λ * step
        ok = c -> _valid(c, dims)
    else
        J = _jacobian(θ, dims, logv)
        φ = _to_phi(θ, dims, logv)
        step = _modify_ai(J' * ev.AI * J) \ (J' * ev.grad)
        big = maximum(abs, step)
        big > maxstep && (step .*= maxstep / big)
        cand = λ -> _from_phi(_floor_phi!(φ + λ * step, dims, logv), dims, logv)
        ok = c -> true
    end
    λ = 1.0
    for _ in 0:maxhalve
        c = cand(λ)
        ok(c) && _reml_eval(ctx, c; full = false).logL ≥ ev.logL && return c, λ
        λ /= 2
    end
    em = _pars([ev.em_R; ev.em_G])
    K === nothing && return em, 0.0
    em = _project(em, θ, ctx.cons)
    _valid(em, dims) && _reml_eval(ctx, em; full = false).logL ≥ ev.logL ? (em, 0.0) : (θ, 0.0)
end

# ---- assembly --------------------------------------------------------
# The coefficient matrix C keeps one sparsity pattern over all iterates:
# the union of the patterns of W_a′W_b for every residual block and pair
# of its traits, and of every block (a, b) of G₀ₖ⁻¹ ⊗ Kₖ⁻¹. Positions of
# every contribution in that pattern are computed once, so that each
# evaluation only refills values. R⁻¹ (block diagonal by residual group)
# likewise keeps one pattern. The sparse back-end reuses the symbolic
# Cholesky factorisation (CHOLMOD) and extracts C⁻¹ on the pattern of C
# from the selected inverse (SelectedInversion.jl) with a plan made once.

# index into `S.nzval` of every stored entry of `M`, shifted by (ro, co)
function _positions(S::SparseMatrixCSC, M::SparseMatrixCSC, ro::Integer, co::Integer)
    out = Vector{Int}(undef, nnz(M))
    rows = rowvals(S)
    rm = rowvals(M)
    t = 0
    for j in axes(M, 2)
        lo, hi = S.colptr[j+co], S.colptr[j+co+1] - 1
        for p in nzrange(M, j)
            k = searchsortedfirst(rows, rm[p] + ro, lo, hi, Base.Order.Forward)
            (k ≤ hi && rows[k] == rm[p] + ro) || error("entry outside the pattern of C")
            out[t+=1] = k
        end
    end
    out
end

# Σ_t E[pos[t]] v[t]: the trace tr(E M) over the stored entries of M
function _dotpos(E::Vector{Float64}, pos::Vector{Int}, v::Vector{Float64})
    acc = 0.0
    @inbounds @simd for t in eachindex(pos, v)
        acc += E[pos[t]] * v[t]
    end
    acc
end

mutable struct _SparseState
    F::Any          # CHOLMOD factor, reused across evaluations
    plan::Any       # SelectedInversion extraction plan on C[p, p] (supernodal)
    Cp::Any         # pattern of C[p, p], p the factor's permutation
    back::Vector{Int}  # entry t of C[p, p] is entry back[t] of C
end

function _setup(yv, X, W, effects, blocks, zero)
    n = length(yv)
    keepcol = [setdiff(1:size(X, 2), zero); size(X, 2)+1:size(W, 2)]
    Wk = W[:, keepcol]
    neq = size(Wk, 2)
    # M[i] = W_a′W_b + W_b′W_a (W_a′W_a for a = b) of block `ib`, traits
    # a ≥ b (local indices in P), with Mkey[i] = (ib, a, b)
    M = SparseMatrixCSC{Float64,Int}[]
    Mkey = Tuple{Int,Int,Int}[]
    for (ib, b) in enumerate(blocks)
        Ws = [Wk[b.rows[a, :], :] for a in eachindex(b.P)]
        for bb in eachindex(b.P), a in bb:length(b.P)
            Mab = sparse(Ws[a]' * Ws[bb])
            push!(M, a == bb ? Mab : Mab + sparse(Mab'))
            push!(Mkey, (ib, a, bb))
        end
    end
    off = Int[]                                   # first equation of each effect, minus one
    o = neq - sum(size(r.Z, 2) for r in effects; init = 0)
    K = SparseMatrixCSC{Float64,Int}[]
    for r in effects
        push!(off, o)
        push!(K, r.Kinv isa SparseMatrixCSC ? r.Kinv : sparse(Matrix(r.Kinv)))
        o += size(r.Z, 2)
    end
    # pattern of C
    I_, J_ = Int[], Int[]
    for Mi in M
        i, j, _ = findnz(Mi)
        append!(I_, i); append!(J_, j)
    end
    for (k, r) in enumerate(effects)
        q, mk = nlevels(r), ncomponents(r)
        i, j, _ = findnz(K[k])
        for a in 1:mk, b in 1:mk
            append!(I_, i .+ (off[k] + (a - 1) * q)); append!(J_, j .+ (off[k] + (b - 1) * q))
        end
    end
    S = sparse(I_, J_, 1.0, neq, neq)
    posM = [_positions(S, Mi, 0, 0) for Mi in M]
    posK = [[_positions(S, K[k], off[k] + (a - 1) * nlevels(r), off[k] + (b - 1) * nlevels(r))
             for a in 1:ncomponents(r), b in 1:ncomponents(r)] for (k, r) in enumerate(effects)]
    # pattern of R⁻¹; Rslot[t]: index of stored entry t in the
    # concatenation of vec(Q_b), Q_b = R₀[P, P]⁻¹ of block b
    I_, J_, sl = Int[], Int[], Int[]
    qo = 0
    for b in blocks
        np = length(b.P)
        for g in axes(b.rows, 2), bb in 1:np, a in 1:np
            push!(I_, b.rows[a, g]); push!(J_, b.rows[bb, g]); push!(sl, qo + (bb - 1) * np + a)
        end
        qo += np^2
    end
    Rp = sparse(I_, J_, collect(1:length(I_)), n, n)
    Rslot = sl[nonzeros(Rp)]
    (; keepcol, Wk, M, Mkey, off, K, S, posM, posK, Rp, Rslot,
     state = _SparseState(nothing, nothing, nothing, Int[]))
end

# C at the residual inverses Q (per block) and G0
function _fill_C(st, Q, G0)
    S = st.S
    nz = zeros(nnz(S))
    for (i, (ib, a, b)) in enumerate(st.Mkey)
        w = Q[ib][a, b]
        v = nonzeros(st.M[i])
        pos = st.posM[i]
        @inbounds for t in eachindex(pos)
            nz[pos[t]] += w * v[t]
        end
    end
    for k in eachindex(st.K)
        Gi = inv(Symmetric(G0[k]))
        v = nonzeros(st.K[k])
        for a in axes(Gi, 1), b in axes(Gi, 2)
            g = Gi[a, b]
            pos = st.posK[k][a, b]
            @inbounds for t in eachindex(pos)
                nz[pos[t]] += g * v[t]
            end
        end
    end
    SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, nz)
end

function _sparse_factor!(st, C)
    if st.state.F === nothing
        st.state.F = cholesky(Symmetric(C); check = false)
    else
        cholesky!(st.state.F, Symmetric(C); check = false)
    end
    st.state.F
end

# C⁻¹ on the pattern of C (values aligned with S.nzval). The selected
# inverse is taken in the factor's ordering p and mapped back through a
# permutation of entry indices made once; building the extraction plan on
# C[p, p] is much cheaper than SelectedInversion's own depermutation.
function _sparse_cinv(st, F, C)
    Z, p = selinv(F; depermute = false)
    stt = st.state
    if stt.Cp === nothing
        tag = SparseMatrixCSC(C.m, C.n, C.colptr, C.rowval, collect(1:nnz(C)))[p, p]
        stt.back = nonzeros(tag)
        stt.Cp = SparseMatrixCSC(C.m, C.n, tag.colptr, tag.rowval, zeros(nnz(C)))
        Z isa SupernodalMatrix && (stt.plan = selinv_extract_setup(Z, stt.Cp))
    end
    Ep = stt.Cp
    if Z isa SupernodalMatrix
        selinv_extract!(Ep, Z, stt.plan)
    else
        selinv_extract!(Ep, Z, stt.Cp)
    end
    Ev = Vector{Float64}(undef, nnz(C))
    Ev[stt.back] = nonzeros(Ep)
    Ev
end

# dense C⁻¹ gathered on the pattern of C
function _gather(Ci, S)
    rows = rowvals(S)
    Ev = Vector{Float64}(undef, nnz(S))
    for j in axes(S, 2), p in nzrange(S, j)
        Ev[p] = Ci[rows[p], j]
    end
    Ev
end

# ---- evaluation ------------------------------------------------------
# log-likelihood at θ; with `full`, also gradient, AI matrix (covariance
# scale), EM updates and solutions
function _reml_eval(ctx, θ; full::Bool)
    (; yv, W, effects, blocks, nR, dims, logK, hom, p, be, store) = ctx
    mats = _unpack(θ, dims)
    all(M -> isposdef(Symmetric(M)), mats) || return (logL = -Inf,)
    R, G0 = mats[1:nR], mats[nR+1:end]
    n = length(yv)
    Q = [Matrix(inv(Symmetric(R[b.c][b.P, b.P]))) for b in blocks]
    ldR = sum(size(b.rows, 2) * logdet(Symmetric(R[b.c][b.P, b.P])) for b in blocks)
    Rp = store.Rp
    Rinv = SparseMatrixCSC(n, n, Rp.colptr, Rp.rowval, reduce(vcat, vec.(Q))[store.Rslot])
    ry = Rinv * yv
    rhs = store.Wk' * ry
    C = _fill_C(store, Q, G0)
    if be == :sparse
        F = _sparse_factor!(store, C)
        issuccess(F) || return (logL = -Inf,)
    else
        F = cholesky(Symmetric(Matrix(C)); check = false)
        issuccess(F) || return (logL = -Inf,)
    end
    xk = F \ rhs
    x = zeros(size(W, 2))
    x[store.keepcol] = xk
    e = yv - store.Wk * xk
    ldG = sum(size(r.Kinv, 1) * logdet(G0[k]) + size(G0[k], 1) * logK[k]
              for (k, r) in enumerate(effects); init = 0.0)
    logL = -0.5 * (dot(ry, e) + logdet(F) + ldR + ldG)
    full || return (logL = logL,)
    # C⁻¹ on the pattern of C; solves with C
    if be == :sparse
        Ev = _sparse_cinv(store, F, C)
        quadW = RF -> (WRF = store.Wk' * RF; WRF' * (F \ WRF))
    else
        Ci = inv(F)
        Ev = _gather(Ci, store.S)
        quadW = RF -> (WRF = store.Wk' * RF; WRF' * (Ci * WRF))
    end
    # traces: T[ib][a,b] = tr(C⁻¹W_a′W_b) per residual block,
    # Ts[k][a,b] = tr(C^{ab}K_k⁻¹) per random effect
    T = [zeros(length(b.P), length(b.P)) for b in blocks]
    for (i, (ib, a, b)) in enumerate(store.Mkey)
        v = _dotpos(Ev, store.posM[i], nonzeros(store.M[i]))
        T[ib][a, b] = T[ib][b, a] = a == b ? v : v / 2
    end
    Ts = [[_dotpos(Ev, store.posK[k][a, b], nonzeros(store.K[k]))
           for a in 1:ncomponents(r), b in 1:ncomponents(r)] for (k, r) in enumerate(effects)]
    # residual parameters: per block, Σ_g r_g r_g′ (r = R⁻¹ê) and Σ_g ê_g ê_g′
    re = Rinv * e
    Rr = [(v = re[b.rows]; v * v') for b in blocks]
    Ee = [(v = e[b.rows]; v * v') for b in blocks]
    grad = Float64[]
    fs = Vector{Vector{Float64}}()
    em_R = Matrix{Float64}[]
    for c in 1:nR
        t = dims[c]
        for (i, j) in _vech(t)
            D = zeros(t, t)
            D[i, j] = D[j, i] = 1
            val = 0.0
            f = zeros(n)
            for (ib, b) in enumerate(blocks)
                b.c == c || continue
                DP = D[b.P, b.P]
                iszero(DP) && continue
                Qb = Q[ib]
                val += size(b.rows, 2) * tr(Qb * DP) - sum((Qb * DP * Qb) .* T[ib]) - sum(DP .* Rr[ib])
                f[b.rows] = DP * re[b.rows]
            end
            push!(grad, -0.5 * val)
            push!(fs, f)
        end
        if hom
            push!(em_R, fill(dot(e, yv) / (n - p), 1, 1))
        else
            Rc = R[c]
            acc = zeros(t, t)
            N = 0
            for (ib, b) in enumerate(blocks)
                b.c == c || continue
                P = b.P
                ng = size(b.rows, 2)
                B = zeros(t, length(P))
                B[P, :] = I(length(P))
                Mi = setdiff(1:t, P)
                if !isempty(Mi)
                    B[Mi, :] = Rc[Mi, P] / Symmetric(Rc[P, P])
                    acc[Mi, Mi] .+= ng .* (Rc[Mi, Mi] - Rc[Mi, P] * (Symmetric(Rc[P, P]) \ Rc[P, Mi]))
                end
                acc .+= B * (Ee[ib] + T[ib]) * B'
                N += ng
            end
            push!(em_R, (acc + acc') ./ 2N)
        end
    end
    # random effects
    em_G = Matrix{Float64}[]
    blocksU = [store.keepcol[store.off[k]+1:store.off[k]+size(r.Z, 2)] for (k, r) in enumerate(effects)]
    for (k, r) in enumerate(effects)
        q, mk = nlevels(r), ncomponents(r)
        U = reshape(x[blocksU[k]], q, mk)
        Gi = inv(Symmetric(G0[k]))
        Tk = Ts[k]
        S = U' * (r.Kinv * U)                     # Û'K⁻¹Û
        push!(em_G, (S + Tk) ./ q)
        for (i, j) in _vech(mk)
            D = zeros(mk, mk)
            D[i, j] = D[j, i] = 1
            trP = q * tr(Gi * D) - tr(Tk * Gi * D * Gi)
            quad = tr(Gi * D * Gi * S)
            push!(grad, -0.5 * (trP - quad))
            push!(fs, r.Z * vec(U * (D * Gi)'))
        end
    end
    # average information: AI_ij = ½ f_i' P f_j, P f = R⁻¹f − R⁻¹W C⁻¹ W'R⁻¹ f
    Fm = reduce(hcat, fs)
    RF = Rinv * Fm
    AI = 0.5 .* (Fm' * RF .- quadW(RF))
    (logL = logL, grad = grad, AI = Symmetric(AI), em_R = em_R, em_G = em_G, x = x)
end

# log-determinant of a symmetric positive definite matrix; a sparse matrix
# is factored sparsely, so K⁻¹ of a large pedigree is never made dense
_logdet_spd(K::SparseMatrixCSC) = logdet(cholesky(Symmetric(K)))
_logdet_spd(K::AbstractMatrix) = logdet(cholesky(Symmetric(Matrix(K))))
