"""
    ResidualStructure

Residual (co)variance structure for [`reml`](@ref):
[`Homogeneous`](@ref) or [`HetClasses`](@ref).
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

_variances(r::Homogeneous) = [r.σ²]
_variances(r::HetClasses) = copy(r.σ²)
_classes(::Homogeneous, n) = ones(Int, n)
function _classes(r::HetClasses, n)
    length(r.class) == n || throw(DimensionMismatch("$(length(r.class)) classes for $n records"))
    r.class
end

"""
    REMLResult

Output of [`reml`](@ref), all at the returned estimates:

- `σ²e`: the residual variance (homogeneous residuals; `NaN` for
  [`HetClasses`](@ref)), `residual`: the residual variance of every
  class;
- `G0`: the covariance matrices of the random effects;
- `logL`: the restricted log-likelihood (Mrode's form, without
  constants, as in WOMBAT);
- `AIinv`: the inverse of the average-information matrix on the
  covariance scale (approximate sampling covariances); parameter order:
  residual variances, then the lower triangle of each `G0`, by columns;
- `history`: parameters per iterate (one column each, starting with the
  initial values), `logLs`: their log-likelihoods;
- `iterations`, and the mixed-model `solutions`.
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
end

# lower-triangle index pairs of an m × m covariance matrix
_vech(m) = [(i, j) for j in 1:m for i in j:m]

_pars(s, G0) = [s; reduce(vcat, [[G[i, j] for (i, j) in _vech(size(G, 1))] for G in G0]; init = Float64[])]

function _unpack(θ, ns, G0)
    s = θ[1:ns]
    k = ns
    G = map(G0) do g
        mk = size(g, 1)
        new = zeros(mk, mk)
        for (i, j) in _vech(mk)
            k += 1
            new[i, j] = new[j, i] = θ[k]
        end
        new
    end
    s, G
end

function _valid(θ, ns, G0)
    s, G = _unpack(θ, ns, G0)
    all(>(0), s) && all(isposdef, G)
end

# log-Cholesky parameters φ: log σ² for residuals; for G = LL′ the lower
# triangle of L by columns with log-transformed diagonal
function _to_phi(θ, ns, G0)
    s, G = _unpack(θ, ns, G0)
    φ = log.(s)
    for g in G
        L = cholesky(Symmetric(g)).L
        append!(φ, [i == j ? log(L[i, i]) : L[i, j] for (i, j) in _vech(size(g, 1))])
    end
    φ
end

function _from_phi(φ, ns, G0)
    θ = exp.(φ[1:ns])
    k = ns
    for g in G0
        m = size(g, 1)
        L = zeros(m, m)
        for (i, j) in _vech(m)
            k += 1
            L[i, j] = i == j ? exp(φ[k]) : φ[k]
        end
        G = L * L'
        append!(θ, [G[i, j] for (i, j) in _vech(m)])
    end
    θ
end

# keep every Cholesky pivot at least `rel` times the largest one of its
# matrix (WOMBAT's limit for "small pivots"), so that estimates on the
# boundary of the parameter space stay numerically positive definite
function _floor_phi!(φ, ns, G0; rel = 1e-4)
    k = ns
    for g in G0
        pr = _vech(size(g, 1))
        d = [k + c for (c, (i, j)) in enumerate(pr) if i == j]
        lo = maximum(φ[d]) + log(rel)
        φ[d] .= max.(φ[d], lo)
        k += length(pr)
    end
    φ
end

# Jacobian ∂θ/∂φ (block diagonal)
function _jacobian(θ, ns, G0)
    s, G = _unpack(θ, ns, G0)
    blocks = Matrix{Float64}[Matrix(Diagonal(s))]
    for g in G
        m = size(g, 1)
        L = cholesky(Symmetric(g)).L
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
         zero = Int[], tol = 1e-8, ltol = nothing, maxiter = 200, maxhalve = 10,
         maxstep = 5.0) -> REMLResult

Restricted maximum likelihood (Patterson and Thompson, 1971) for the
univariate mixed model

```math
y = Xb + \\sum_k Z_k u_k + e, \\quad \\operatorname{var}(u_k) = G_{0k} \\otimes K_k,
\\quad \\operatorname{var}(e) = R,
```

where each [`RandomEffect`](@ref) in `effects` supplies ``Z_k``,
``K_k^{-1}`` and the starting ``G_{0k}`` (so correlated direct–maternal
effects or random-regression covariance matrices are estimated too).
The residual structure is given either as `σ²e` (starting value of a
homogeneous variance) or as `residual`, a [`Homogeneous`](@ref) or
[`HetClasses`](@ref) with starting values. `zero` lists constrained
fixed effects (as in [`solve_mme`](@ref)).

`method`:
- `:ai` – average information (Gilmour et al., 1995): first derivatives
  from MME solutions and traces of ``C^{-1}``, and
  ``AI_{ij} = \\tfrac12 f_i'Pf_j`` from working variates
  ``f_i = (\\partial V/\\partial\\theta_i)Py`` (``D_c e/\\sigma^2_c`` for
  residual class ``c``, ``Z_k(\\partial G_{0k} G_{0k}^{-1} \\otimes
  I)\\hat u_k`` for ``G_{0k}``). The AI matrix is made safely positive
  definite by adding ``\\tau I``, ``\\tau = \\max(0, \\delta - \\lambda_{min})``,
  ``\\delta = 3 \\times 10^{-6}\\lambda_{max}`` (WOMBAT's default
  modification). A step must not decrease the log-likelihood; otherwise
  it is halved up to `maxhalve` times, then replaced by an EM step.
- `:em` – expectation maximization: for homogeneous residuals
  ``(n - p)\\sigma^2_e = e'y``, for classes ``n_c\\sigma^2_c = e_c'e_c +
  \\operatorname{tr}(C^{-1}W_c'W_c)``, and ``G_{0k} = (\\hat
  U_k'K_k^{-1}\\hat U_k + T_k)/q_k`` with ``T_k[i,j] =
  \\operatorname{tr}(C^{k_ik_j}K_k^{-1})``.

`param` (AI only): `:natural` takes AI steps on the covariances
themselves (steps leaving the parameter space are halved);
`:cholesky` on log residual variances and the Cholesky factors of the
``G_{0k}`` with log-transformed diagonals, which keeps every iterate
positive definite. The average information is computed on the
covariance scale and transformed with the Jacobian (as in WOMBAT).
With `:cholesky`, no element of a step may exceed `maxstep` (the step is
scaled down), and every Cholesky pivot is kept at least ``10^{-4}``
times the largest pivot of its matrix, so that estimates on the boundary
of the parameter space remain numerically positive definite.

Iteration stops when the largest change in a covariance parameter is
below `tol` times the largest parameter (at least 1); or, if `ltol` is
given, when an AI step increases the log-likelihood by less than `ltol`
(WOMBAT uses ``5 \\times 10^{-4}``); or after `maxiter` iterates.
Near the boundary of the parameter space the likelihood can be very
flat in some directions, and the parameter criterion alone may need
many iterates.

The log-likelihood is ``-\\tfrac12[y'Py + \\log|V| + \\log|X'V^{-1}X|]``
computed through ``\\log|C| + \\log|R| + \\sum_k \\log|G_k|``. ``C^{-1}``
is formed densely, which limits the model to a few thousand equations.
"""
function reml(y::AbstractVector, X::AbstractMatrix, effects::AbstractVector{<:RandomEffect};
              σ²e = nothing, residual = nothing, method::Symbol = :ai, param::Symbol = :natural,
              zero = Int[], tol = 1e-8, ltol = nothing, maxiter = 200, maxhalve = 10,
              maxstep = 5.0)
    method in (:ai, :em) || throw(ArgumentError("method must be :ai or :em"))
    param in (:natural, :cholesky) || throw(ArgumentError("param must be :natural or :cholesky"))
    (σ²e === nothing) == (residual === nothing) &&
        throw(ArgumentError("give exactly one of σ²e or residual"))
    res = residual === nothing ? Homogeneous(σ²e) : residual
    hom = res isa Homogeneous
    n = length(y)
    yv = Vector{Float64}(y)
    cls = _classes(res, n)
    s = _variances(res)
    ns = length(s)
    W = SparseMatrixCSC{Float64,Int}(hcat(sparse(X), (r.Z for r in effects)...))
    G0 = [copy(r.G0) for r in effects]
    # log|K_k| from K_k⁻¹ (constant over iterations)
    logK = [-logdet(cholesky(Symmetric(Matrix(r.Kinv)))) for r in effects]
    # W_c′W_c per residual class (traces for the residual derivatives)
    WW = [sparse(W[cls.==c, :]' * W[cls.==c, :]) for c in 1:ns]
    nc = [count(==(c), cls) for c in 1:ns]
    ctx = (; yv, X, W, effects, cls, ns, nc, WW, logK, zero, hom)

    θ = _pars(s, G0)
    ev = _reml_eval(ctx, θ, G0; full = true)
    hist = [θ]
    logLs = [ev.logL]
    it = 0
    while it < maxiter
        θn, ai = method == :em ? (_pars(ev.em_s, ev.em_G), false) :
                 _ai_step(ctx, θ, G0, ev, param, maxhalve, maxstep)
        evn = _reml_eval(ctx, θn, G0; full = true)
        conv = maximum(abs, θn - θ) < tol * max(1, maximum(abs, θ)) ||
               (ltol !== nothing && ai && evn.logL - ev.logL < ltol)
        θ, ev = θn, evn
        push!(hist, θ)
        push!(logLs, ev.logL)
        it += 1
        conv && break
    end
    s, G = _unpack(θ, ns, G0)
    REMLResult(hom ? s[1] : NaN, G, ev.logL, inv(ev.AI), reduce(hcat, hist), logLs, it, ev.x, s)
end

# WOMBAT's default modification of the AI matrix: add τI so that the
# smallest eigenvalue is at least 3e-6 times the largest
function _modify_ai(A)
    λ = eigvals(Symmetric(Matrix(A)))
    τ = max(0.0, 3e-6 * λ[end] - λ[1])
    Symmetric(Matrix(A) + τ * I)
end

# AI step from θ, halved until the log-likelihood does not decrease; EM
# step if no step length does. Returns the new parameters and whether an
# AI step was taken.
function _ai_step(ctx, θ, G0, ev, param, maxhalve, maxstep)
    if param == :natural
        step = _modify_ai(ev.AI) \ ev.grad
        cand = λ -> θ + λ * step
        ok = c -> _valid(c, ctx.ns, G0)
    else
        J = _jacobian(θ, ctx.ns, G0)
        φ = _to_phi(θ, ctx.ns, G0)
        step = _modify_ai(J' * ev.AI * J) \ (J' * ev.grad)
        big = maximum(abs, step)
        big > maxstep && (step .*= maxstep / big)
        cand = λ -> _from_phi(_floor_phi!(φ + λ * step, ctx.ns, G0), ctx.ns, G0)
        ok = c -> true
    end
    λ = 1.0
    for _ in 0:maxhalve
        c = cand(λ)
        ok(c) && _reml_eval(ctx, c, G0; full = false).logL ≥ ev.logL && return c, true
        λ /= 2
    end
    _pars(ev.em_s, ev.em_G), false
end

# log-likelihood at θ; with `full`, also gradient, AI matrix (covariance
# scale), EM updates and solutions
function _reml_eval(ctx, θ, G0ref; full::Bool)
    (; yv, X, W, effects, cls, ns, nc, WW, logK, zero, hom) = ctx
    s, G0 = _unpack(θ, ns, G0ref)
    (all(>(0), s) && all(isposdef, G0)) || return (logL = -Inf,)
    n = length(yv)
    rinv = 1 ./ s[cls]
    effs = [RandomEffect(r.Z, r.Kinv, G0[k], r.name) for (k, r) in enumerate(effects)]
    m = mme(yv, X, effs; Rinv = rinv)
    nt = length(m.rhs)
    keep = _keep(nt, zero)
    p = count(≤(size(X, 2)), keep)
    F = cholesky(Symmetric(Matrix(m.lhs[keep, keep])); check = false)
    issuccess(F) || return (logL = -Inf,)
    x = zeros(nt)
    x[keep] = F \ m.rhs[keep]
    e = yv - W * x
    ldG = sum(size(r.Kinv, 1) * logdet(G0[k]) + size(G0[k], 1) * logK[k]
              for (k, r) in enumerate(effects); init = 0.0)
    ldR = sum(nc[c] * log(s[c]) for c in 1:ns)
    logL = -0.5 * (dot(rinv .* yv, e) + logdet(F) + ldR + ldG)
    full || return (logL = logL,)
    Ci = zeros(nt, nt)
    Ci[keep, keep] = inv(F)
    # residual parameters
    grad = Float64[]
    fs = Vector{Vector{Float64}}()
    em_s = zeros(ns)
    for c in 1:ns
        I_, J_, V_ = findnz(WW[c])
        trc = sum(Ci[i, j] * v for (i, j, v) in zip(I_, J_, V_); init = 0.0)
        ec = e .* (cls .== c)
        push!(grad, -0.5 * (nc[c] / s[c] - trc / s[c]^2 - dot(ec, ec) / s[c]^2))
        push!(fs, ec ./ s[c])
        em_s[c] = hom ? dot(e, yv) / (n - p) : (dot(ec, ec) + trc) / nc[c]
    end
    # random effects
    em_G = Matrix{Float64}[]
    for (k, r) in enumerate(effects)
        rng = m.blocks[k+1]
        q, mk = nlevels(r), ncomponents(r)
        U = reshape(x[rng], q, mk)
        Gi = inv(Symmetric(G0[k]))
        Ck_ = Ci[rng, rng]
        T = [tr(Ck_[(a-1)*q+1:a*q, (b-1)*q+1:b*q] * r.Kinv) for a in 1:mk, b in 1:mk]
        S = U' * (r.Kinv * U)                     # Û'K⁻¹Û
        push!(em_G, (S + T) ./ q)
        for (i, j) in _vech(mk)
            D = zeros(mk, mk)
            D[i, j] = D[j, i] = 1
            trP = q * tr(Gi * D) - tr(T * Gi * D * Gi)
            quad = tr(Gi * D * Gi * S)
            push!(grad, -0.5 * (trP - quad))
            push!(fs, r.Z * vec(U * (D * Gi)'))
        end
    end
    # average information: AI_ij = ½ f_i' P f_j, P f = R⁻¹f − R⁻¹W C⁻¹ W'R⁻¹ f
    Fm = reduce(hcat, fs)
    RF = rinv .* Fm
    WRF = W' * RF
    AI = 0.5 .* (Fm' * RF .- WRF' * (Ci * WRF))
    (logL = logL, grad = grad, AI = Symmetric(AI), em_s = em_s, em_G = em_G, x = x)
end
