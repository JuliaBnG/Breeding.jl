"""
    REMLResult

Output of [`reml`](@ref): the residual variance `σ²e`, the estimated
covariance matrices `G0` of the random effects, the restricted
log-likelihood `logL` (Mrode's form, without constants), the inverse of
the average-information matrix `AIinv` at the last iterate (approximate
sampling covariances; parameter order `σ²e`, then the lower triangle of
each `G0`, column-major), the history of parameters and log-likelihoods
(one column per iterate, starting with the initial values), and the
final mixed-model solutions.
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
end

# lower-triangle index pairs of an m × m covariance matrix
_vech(m) = [(i, j) for j in 1:m for i in j:m]

"""
    reml(y, X, effects; σ²e, method = :ai, zero = Int[], tol = 1e-8, maxiter = 200)
        -> REMLResult

Restricted maximum likelihood (Patterson and Thompson, 1971) for the
univariate mixed model

```math
y = Xb + \\sum_k Z_k u_k + e, \\quad \\operatorname{var}(u_k) = G_{0k} \\otimes K_k,
\\quad \\operatorname{var}(e) = I\\sigma^2_e,
```

where each [`RandomEffect`](@ref) in `effects` supplies ``Z_k``,
``K_k^{-1}`` and the starting ``G_{0k}`` (so correlated direct–maternal
effects or random-regression covariance matrices are estimated too).
`σ²e` is the starting residual variance and `zero` lists constrained
fixed effects (as in [`solve_mme`](@ref)).

`method`:
- `:ai` – average information (Gilmour et al., 1995): first derivatives
  from MME solutions and traces of ``C^{-1}``, and
  ``AI_{ij} = \\tfrac12 f_i'Pf_j`` from working variates
  ``f_i = (\\partial V/\\partial\\theta_i)Py`` (``e/\\sigma^2_e`` for the
  residual, ``Z_k(\\partial G_{0k} G_{0k}^{-1} \\otimes I)\\hat u_k`` for
  ``G_{0k}``); steps are halved if a matrix leaves the positive
  definite cone;
- `:em` – expectation maximization: ``(n - p)\\sigma^2_e = e'y`` and
  ``G_{0k} = (\\hat U_k'K_k^{-1}\\hat U_k + T_k)/q_k`` with
  ``T_k[i,j] = \\operatorname{tr}(C^{k_ik_j}K_k^{-1})``.

The log-likelihood is ``-\\tfrac12[y'Py + \\log|V| + \\log|X'V^{-1}X|]``
computed through ``\\log|C| + \\log|R| + \\sum_k \\log|G_k|``. ``C^{-1}`` is
formed densely, which limits the model to a few thousand equations.
"""
function reml(y::AbstractVector, X::AbstractMatrix, effects::AbstractVector{<:RandomEffect};
              σ²e::Real, method::Symbol = :ai, zero = Int[], tol = 1e-8, maxiter = 200)
    method in (:ai, :em) || throw(ArgumentError("method must be :ai or :em"))
    n = length(y)
    yv = Vector{Float64}(y)
    W = SparseMatrixCSC{Float64,Int}(hcat(sparse(X), (r.Z for r in effects)...))
    G0 = [copy(r.G0) for r in effects]
    ve = Float64(σ²e)
    # log|K_k| from K_k⁻¹ (constant over iterations)
    logK = [-logdet(cholesky(Symmetric(Matrix(r.Kinv)))) for r in effects]
    pars(ve, G0) = [ve; reduce(vcat, [[G[i, j] for (i, j) in _vech(size(G, 1))] for G in G0])]
    hist = [pars(ve, G0)]
    logLs = Float64[]
    AIinv = zeros(0, 0)
    x = Float64[]
    it = 0
    while true
        effs = [RandomEffect(r.Z, r.Kinv, G0[k], r.name) for (k, r) in enumerate(effects)]
        m = mme(yv, X, effs; Rinv = fill(1 / ve, n))
        nt = length(m.rhs)
        keep = _keep(nt, zero)
        p = count(≤(size(X, 2)), keep)
        Ck = Symmetric(Matrix(m.lhs[keep, keep]))
        F = cholesky(Ck)
        Ci = zeros(nt, nt)
        Ci[keep, keep] = inv(F)
        x = zeros(nt)
        x[keep] = F \ m.rhs[keep]
        e = yv - W * x
        # log-likelihood
        ldG = sum(size(r.Kinv, 1) * logdet(G0[k]) + size(G0[k], 1) * logK[k]
                  for (k, r) in enumerate(effects); init = 0.0)
        push!(logLs, -0.5 * (dot(yv, e) / ve + logdet(F) + n * log(ve) + ldG))
        # gradient, working variates and EM updates
        grad = Float64[]
        fs = Vector{Vector{Float64}}()
        dfu = 0.0                                     # Σ (dim_k − tr(C^{kk} G_k⁻¹))
        emG = Matrix{Float64}[]
        for (k, r) in enumerate(effects)
            rng = m.blocks[k+1]
            q, mk = nlevels(r), ncomponents(r)
            U = reshape(x[rng], q, mk)
            Gi = inv(Symmetric(G0[k]))
            Ck_ = Ci[rng, rng]
            T = [tr(Ck_[(a-1)*q+1:a*q, (b-1)*q+1:b*q] * r.Kinv) for a in 1:mk, b in 1:mk]
            S = U' * (r.Kinv * U)                     # Û'K⁻¹Û
            dfu += q * mk - tr(T * Gi)
            push!(emG, (S + T) ./ q)
            for (i, j) in _vech(mk)
                D = zeros(mk, mk)
                D[i, j] = D[j, i] = 1
                trP = q * tr(Gi * D) - tr(T * Gi * D * Gi)
                quad = tr(Gi * D * Gi * S)
                push!(grad, -0.5 * (trP - quad))
                push!(fs, r.Z * vec(U * (D * Gi)'))
            end
        end
        ge = -0.5 * ((n - p - dfu) / ve - dot(e, e) / ve^2)
        pushfirst!(grad, ge)
        pushfirst!(fs, e ./ ve)
        # average information: AI_ij = ½ f_i' P f_j, P f = R⁻¹f − R⁻¹W C⁻¹ W'R⁻¹ f
        Fm = reduce(hcat, fs)
        WRF = W' * Fm ./ ve
        Sol = Ci * WRF
        AI = 0.5 .* (Fm' * Fm ./ ve .- WRF' * Sol)
        AIinv = inv(Symmetric(AI))
        θ = pars(ve, G0)
        θn = if method == :ai
            step = AIinv * grad
            λ = 1.0
            while true
                cand = θ + λ * step
                _valid(cand, G0) && break
                λ /= 2
                λ < 1e-6 && (λ = 0.0; break)
            end
            θ + λ * step
        else
            pars(dot(e, yv) / (n - p), emG)
        end
        ve, G0 = _unpack(θn, G0)
        push!(hist, θn)
        it += 1
        (maximum(abs, θn - θ) < tol * max(1, maximum(abs, θ)) || it ≥ maxiter) && break
    end
    REMLResult(ve, G0, logLs[end], AIinv, reduce(hcat, hist), logLs, it, x)
end

function _unpack(θ, G0)
    ve = θ[1]
    k = 1
    G = map(G0) do g
        mk = size(g, 1)
        new = zeros(mk, mk)
        for (i, j) in _vech(mk)
            k += 1
            new[i, j] = new[j, i] = θ[k]
        end
        new
    end
    ve, G
end

_valid(θ, G0) = θ[1] > 0 && all(isposdef, last(_unpack(θ, G0)))
