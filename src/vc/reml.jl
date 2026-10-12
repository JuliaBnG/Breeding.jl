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
         zero = Int[], tol = 1e-8, ltol = nothing, ntol = nothing, maxiter = 200,
         maxhalve = 10, maxstep = 5.0, backend = :auto) -> REMLResult

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

`backend`: `:dense` forms ``C^{-1}`` densely (small models); `:sparse`
keeps one sparsity pattern of ``C`` over all iterates, reuses the
symbolic CHOLMOD factorisation, takes the traces from the selected
inverse of ``C`` on its own pattern (Takahashi equations, via
SelectedInversion.jl) and the AI matrix from solves with the factor;
`:auto` (default) uses `:sparse` above 2000 equations. Dense
relationship matrices (e.g. ``G^{-1}``) make the corresponding blocks
of ``C`` and its factor dense.
"""
function reml(y::AbstractVector, X::AbstractMatrix, effects::AbstractVector{<:RandomEffect};
              σ²e = nothing, residual = nothing, method::Symbol = :ai, param::Symbol = :natural,
              zero = Int[], tol = 1e-8, ltol = nothing, ntol = nothing, maxiter = 200,
              maxhalve = 10, maxstep = 5.0, backend::Symbol = :auto)
    method in (:ai, :em) || throw(ArgumentError("method must be :ai or :em"))
    param in (:natural, :cholesky) || throw(ArgumentError("param must be :natural or :cholesky"))
    backend in (:auto, :dense, :sparse) || throw(ArgumentError("backend must be :auto, :dense or :sparse"))
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
    nc = [count(==(c), cls) for c in 1:ns]
    neq = size(W, 2) - length(zero)
    be = backend == :auto ? (neq > 2000 ? :sparse : :dense) : backend
    store = be == :sparse ? _sparse_setup(yv, X, W, effects, cls, ns, zero) :
            (WW = [sparse(W[cls.==c, :]' * W[cls.==c, :]) for c in 1:ns],)
    ctx = (; yv, X, W, effects, cls, ns, nc, logK, zero, hom, be, store)

    θ = _pars(s, G0)
    ev = _reml_eval(ctx, θ, G0; full = true)
    hist = [θ]
    logLs = [ev.logL]
    steps = Float64[]
    it = 0
    while it < maxiter
        θn, λ = method == :em ? (_pars(ev.em_s, ev.em_G), 0.0) :
                _ai_step(ctx, θ, G0, ev, param, maxhalve, maxstep)
        evn = _reml_eval(ctx, θn, G0; full = true)
        conv = maximum(abs, θn - θ) < tol * max(1, maximum(abs, θ)) ||
               (ltol !== nothing && λ > 0 && evn.logL - ev.logL < ltol &&
                (ntol === nothing || _gain(evn, θn, ns, G0, param) < ntol))
        θ, ev = θn, evn
        push!(hist, θ)
        push!(logLs, ev.logL)
        push!(steps, λ)
        it += 1
        conv && break
    end
    s, G = _unpack(θ, ns, G0)
    REMLResult(hom ? s[1] : NaN, G, ev.logL, inv(ev.AI), reduce(hcat, hist), logLs, it, ev.x, s, steps)
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
# increase in log L expected from a (modified) AI step: ½ g′H⁻¹g
function _gain(ev, θ, ns, G0, param)
    g, H = ev.grad, ev.AI
    if param == :cholesky
        J = _jacobian(θ, ns, G0)
        g, H = J' * g, J' * H * J
    end
    dot(g, _modify_ai(H) \ g) / 2
end

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
        ok(c) && _reml_eval(ctx, c, G0; full = false).logL ≥ ev.logL && return c, λ
        λ /= 2
    end
    _pars(ev.em_s, ev.em_G), 0.0
end

# ---- sparse back-end -------------------------------------------------
# The coefficient matrix C keeps one sparsity pattern over all iterates:
# the union of the patterns of W_c′W_c (residual classes) and of every
# block (a, b) of G₀ₖ⁻¹ ⊗ Kₖ⁻¹. Positions of every contribution in that
# pattern are computed once, so that each evaluation only refills values,
# the symbolic Cholesky factorisation (CHOLMOD) is reused, and C⁻¹ on the
# pattern of C is extracted from the selected inverse (SelectedInversion.jl)
# with a plan made once.

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

function _sparse_setup(yv, X, W, effects, cls, ns, zero)
    keepcol = [setdiff(1:size(X, 2), zero); size(X, 2)+1:size(W, 2)]
    Wk = W[:, keepcol]
    neq = size(Wk, 2)
    WW = SparseMatrixCSC{Float64,Int}[]
    rhs = Vector{Float64}[]
    for c in 1:ns
        rows = findall(==(c), cls)
        Wc = Wk[rows, :]
        push!(WW, sparse(Wc' * Wc))
        push!(rhs, Wc' * yv[rows])
    end
    off = Int[]                                   # first equation of each effect, minus one
    o = length(keepcol) - sum(size(r.Z, 2) for r in effects; init = 0)
    K = SparseMatrixCSC{Float64,Int}[]
    for r in effects
        push!(off, o)
        push!(K, r.Kinv isa SparseMatrixCSC ? r.Kinv : sparse(Matrix(r.Kinv)))
        o += size(r.Z, 2)
    end
    # pattern of C
    I_, J_ = Int[], Int[]
    for M in WW
        i, j, _ = findnz(M)
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
    I_ = J_ = nothing
    posWW = [_positions(S, M, 0, 0) for M in WW]
    posK = [[_positions(S, K[k], off[k] + (a - 1) * nlevels(r), off[k] + (b - 1) * nlevels(r))
             for a in 1:ncomponents(r), b in 1:ncomponents(r)] for (k, r) in enumerate(effects)]
    (; keepcol, Wk, WW, rhs, off, K, S, posWW, posK, state = _SparseState(nothing, nothing, nothing, Int[]))
end

# factorise C at (s, G0); returns (factor or nothing, solutions in Wk coordinates)
function _sparse_factor(st, s, G0)
    S = st.S
    nz = zeros(nnz(S))
    for c in eachindex(st.WW)
        w = 1 / s[c]
        v = nonzeros(st.WW[c])
        pos = st.posWW[c]
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
    C = SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, nz)
    if st.state.F === nothing
        st.state.F = cholesky(Symmetric(C); check = false)
    else
        cholesky!(st.state.F, Symmetric(C); check = false)
    end
    F = st.state.F
    issuccess(F) || return nothing, Float64[], C
    rhs = sum(st.rhs[c] ./ s[c] for c in eachindex(st.rhs))
    F, F \ rhs, C
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

# log-likelihood at θ; with `full`, also gradient, AI matrix (covariance
# scale), EM updates and solutions
function _reml_eval(ctx, θ, G0ref; full::Bool)
    (; yv, X, W, effects, cls, ns, nc, logK, zero, hom, be, store) = ctx
    s, G0 = _unpack(θ, ns, G0ref)
    (all(>(0), s) && all(isposdef, G0)) || return (logL = -Inf,)
    n = length(yv)
    rinv = 1 ./ s[cls]
    nt = size(W, 2)
    x = zeros(nt)
    if be == :sparse
        F, xk, C = _sparse_factor(store, s, G0)
        F === nothing && return (logL = -Inf,)
        x[store.keepcol] = xk
        logdetC = logdet(F)
        p = length(store.keepcol) - sum(size(r.Z, 2) for r in effects; init = 0)
        blocks = [store.off[k]+1:store.off[k]+size(r.Z, 2) for (k, r) in enumerate(effects)]
        blocks = [store.keepcol[rng] for rng in blocks]
    else
        effs = [RandomEffect(r.Z, r.Kinv, G0[k], r.name) for (k, r) in enumerate(effects)]
        m = mme(yv, X, effs; Rinv = rinv)
        keep = _keep(nt, zero)
        p = count(≤(size(X, 2)), keep)
        Fd = cholesky(Symmetric(Matrix(m.lhs[keep, keep])); check = false)
        issuccess(Fd) || return (logL = -Inf,)
        x[keep] = Fd \ m.rhs[keep]
        logdetC = logdet(Fd)
        blocks = m.blocks[2:end]
    end
    e = yv - W * x
    ldG = sum(size(r.Kinv, 1) * logdet(G0[k]) + size(G0[k], 1) * logK[k]
              for (k, r) in enumerate(effects); init = 0.0)
    ldR = sum(nc[c] * log(s[c]) for c in 1:ns)
    logL = -0.5 * (dot(rinv .* yv, e) + logdetC + ldR + ldG)
    full || return (logL = logL,)
    # traces: trc[c] = tr(C⁻¹W_c′W_c), T[k][a,b] = tr(C^{ab}K_k⁻¹); solves with C
    if be == :sparse
        Ev = _sparse_cinv(store, F, C)
        trcs = [_dotpos(Ev, store.posWW[c], nonzeros(store.WW[c])) for c in 1:ns]
        Ts = [[_dotpos(Ev, store.posK[k][a, b], nonzeros(store.K[k]))
               for a in 1:ncomponents(r), b in 1:ncomponents(r)] for (k, r) in enumerate(effects)]
        quadW = RF -> (WRF = store.Wk' * RF; WRF' * (F \ WRF))
    else
        Ci = zeros(nt, nt)
        Ci[keep, keep] = inv(Fd)
        trcs = map(1:ns) do c
            I_, J_, V_ = findnz(store.WW[c])
            sum(Ci[i, j] * v for (i, j, v) in zip(I_, J_, V_); init = 0.0)
        end
        Ts = map(enumerate(effects)) do (k, r)
            q, mk = nlevels(r), ncomponents(r)
            Ck_ = Ci[blocks[k], blocks[k]]
            [tr(Ck_[(a-1)*q+1:a*q, (b-1)*q+1:b*q] * r.Kinv) for a in 1:mk, b in 1:mk]
        end
        quadW = RF -> (WRF = W' * RF; WRF' * (Ci * WRF))
    end
    # residual parameters
    grad = Float64[]
    fs = Vector{Vector{Float64}}()
    em_s = zeros(ns)
    for c in 1:ns
        trc = trcs[c]
        ec = e .* (cls .== c)
        push!(grad, -0.5 * (nc[c] / s[c] - trc / s[c]^2 - dot(ec, ec) / s[c]^2))
        push!(fs, ec ./ s[c])
        em_s[c] = hom ? dot(e, yv) / (n - p) : (dot(ec, ec) + trc) / nc[c]
    end
    # random effects
    em_G = Matrix{Float64}[]
    for (k, r) in enumerate(effects)
        q, mk = nlevels(r), ncomponents(r)
        U = reshape(x[blocks[k]], q, mk)
        Gi = inv(Symmetric(G0[k]))
        T = Ts[k]
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
    AI = 0.5 .* (Fm' * RF .- quadW(RF))
    (logL = logL, grad = grad, AI = Symmetric(AI), em_s = em_s, em_G = em_G, x = x)
end
