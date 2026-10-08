"""
    ThresholdResult

Solutions of [`threshold_model`](@ref): thresholds `t`, fixed effects
`b`, random effects `u`, the final coefficient matrix `lhs` (its
generalized inverse gives approximate standard errors), the number of
iterations and the history of solutions (one column per iteration).
"""
struct ThresholdResult
    t::Vector{Float64}
    b::Vector{Float64}
    u::Vector{Float64}
    lhs::Matrix{Float64}
    iterations::Int
    history::Matrix{Float64}
end

_ϕ(x) = exp(-x^2 / 2) / sqrt(2π)
_Φ(x) = cdf(Normal(), x)

"""
    threshold_model(N, X, Z, Kinv, λ; t0 = nothing, zero = Int[], tol = 1e-9,
                    maxiter = 100) -> ThresholdResult

Mode of the posterior of thresholds, fixed and random effects for an
ordered categorical trait with a probit link (Gianola and Foulley,
1983; Mrode Section 15.2). `N` is the `s × m` contingency table of
counts (rows: subclasses, columns: ordered categories), `X` and `Z`
relate subclasses to fixed and random effects, and the random effects
have ``\\operatorname{var}(u) = K\\sigma^2_u`` with ``\\lambda =
1/\\sigma^2_u`` (residual variance of the liability is 1).

Each Newton–Raphson (Fisher scoring) round solves Eqn 15.4,

```math
\\begin{bmatrix} Q & L'X & L'Z \\\\ X'L & X'WX & X'WZ \\\\ Z'L & Z'WX & Z'WZ + K^{-1}\\lambda \\end{bmatrix}
\\begin{bmatrix} \\Delta t \\\\ \\Delta b \\\\ \\Delta u \\end{bmatrix} =
\\begin{bmatrix} p \\\\ X'v \\\\ Z'v - K^{-1}\\lambda u \\end{bmatrix},
```

with `Q, L, W, v, p` from Eqns 15.7–15.12. `zero` lists fixed-effect
columns constrained to zero. Starting thresholds default to the probit
of the cumulative proportions; `b = u = 0`. Iteration stops when the
largest increment is below `tol`.
"""
function threshold_model(N::AbstractMatrix, X::AbstractMatrix, Z::AbstractMatrix,
                         Kinv::AbstractMatrix, λ::Real; t0 = nothing, zero = Int[],
                         tol = 1e-9, maxiter = 100)
    s, m = size(N)
    nt = m - 1
    nb, nu = size(X, 2), size(Z, 2)
    n = vec(sum(N; dims = 2))
    Xd, Zd = Matrix{Float64}(X), Matrix{Float64}(Z)
    K = Matrix{Float64}(Kinv)
    t = t0 === nothing ?
        [quantile(Normal(), sum(N[:, 1:k]) / sum(N)) for k in 1:nt] : Vector{Float64}(t0)
    b, u = zeros(nb), zeros(nu)
    keep = setdiff(1:nt+nb+nu, nt .+ zero)
    hist = Vector{Float64}[]
    ϕ = zeros(s, m + 1)              # ϕ[j, k+1] = φ(t_k − a_j), ϕ at t₀ = t_m = 0
    P = zeros(s, m)
    lhs = zeros(nt + nb + nu, nt + nb + nu)
    it = 0
    while it < maxiter
        it += 1
        a = Xd * b + Zd * u
        for j in 1:s
            Fprev = 0.0
            for k in 1:nt
                d = t[k] - a[j]
                ϕ[j, k+1] = _ϕ(d)
                F = _Φ(d)
                P[j, k] = F - Fprev
                Fprev = F
            end
            P[j, m] = 1 - Fprev
        end
        v = [sum(N[j, k] * (ϕ[j, k] - ϕ[j, k+1]) / P[j, k] for k in 1:m) for j in 1:s]
        w = [n[j] * sum((ϕ[j, k] - ϕ[j, k+1])^2 / P[j, k] for k in 1:m) for j in 1:s]
        L = [-n[j] * ϕ[j, k+1] * ((ϕ[j, k+1] - ϕ[j, k]) / P[j, k] -
                                  (ϕ[j, k+2] - ϕ[j, k+1]) / P[j, k+1]) for j in 1:s, k in 1:nt]
        Q = zeros(nt, nt)
        for k in 1:nt
            Q[k, k] = sum(n[j] * ϕ[j, k+1]^2 * (P[j, k] + P[j, k+1]) / (P[j, k] * P[j, k+1])
                          for j in 1:s)
            k < nt && (Q[k+1, k] = Q[k, k+1] =
                           -sum(n[j] * ϕ[j, k+2] * ϕ[j, k+1] / P[j, k+1] for j in 1:s))
        end
        p = [sum((N[j, k] / P[j, k] - N[j, k+1] / P[j, k+1]) * ϕ[j, k+1] for j in 1:s)
             for k in 1:nt]
        WX, WZ = w .* Xd, w .* Zd
        lhs = [Q L'Xd L'Zd; Xd'L Xd'WX Xd'WZ; Zd'L Zd'WX Zd'WZ+λ*K]
        rhs = [p; Xd'v; Zd'v - λ * (K * u)]
        Δ = zeros(length(rhs))
        Δ[keep] = Symmetric(lhs[keep, keep]) \ rhs[keep]
        t .+= Δ[1:nt]
        b .+= Δ[nt+1:nt+nb]
        u .+= Δ[nt+nb+1:end]
        push!(hist, [t; b; u])
        maximum(abs, Δ) < tol && break
    end
    ThresholdResult(t, b, u, lhs, it, reduce(hcat, hist))
end

"""
    category_probabilities(t, η) -> Matrix{Float64}

Probabilities of the ordered categories for linear predictors `η`
(liability means) and thresholds `t`: ``P_k = \\Phi(t_k - \\eta) -
\\Phi(t_{k-1} - \\eta)``. Returns `length(η) × (length(t) + 1)`.
"""
function category_probabilities(t::AbstractVector, η::AbstractVector)
    F = [_Φ(tk - e) for e in η, tk in t]
    [F[:, 1] diff(F; dims = 2) 1 .- F[:, end]]
end

"""
    JointBinaryResult

Solutions of [`joint_binary_quantitative`](@ref): fixed effects `b1` and
random effects `u1` of the quantitative trait, `τ` and `ν` of the
liability adjusted for the residual regression `β` (so that
``u_2 = \\nu + \\beta u_1``), iterations and history.
"""
struct JointBinaryResult
    b1::Vector{Float64}
    u1::Vector{Float64}
    τ::Vector{Float64}
    ν::Vector{Float64}
    β::Float64
    iterations::Int
    history::Matrix{Float64}
end

"""
    joint_binary_quantitative(y1, d, X1, Z1, X2, Z2, Kinv, G, σ²e1, r12;
                              ycentre = mean(y1), zero1 = Int[], zero2 = Int[],
                              tol = 1e-8, maxiter = 50) -> JointBinaryResult

Joint analysis of a quantitative trait (linear model) and a binary
trait (probit threshold model) (Foulley et al., 1983; Mrode Section
15.3). `y1` are the quantitative records, `d` the binary outcomes (`1` =
the category whose probability is ``\\Phi(\\mu)``, e.g. difficult
calving), one record per row of `X1, Z1, X2, Z2`. Random effects have
``\\operatorname{var}(u_1, u_2) = G \\otimes K``; the residual variance of
the liability is 1, ``\\sigma^2_{e1}`` that of trait 1, ``r_{12}`` the
residual correlation.

The liability is adjusted for the residual regression on trait 1,
``\\beta = r_{12}/\\sigma_{e1}/\\sqrt{1 - r_{12}^2}``, giving
``\\tau = b_2 - \\beta b_1``, ``\\nu = u_2 - \\beta u_1`` and
``G_c = [1\\ 0; -\\beta\\ 1]\\,G\\,[1\\ -\\beta; 0\\ 1]``; the liability mean
of record `j` is ``\\mu_j = x_j'\\tau + z_j'\\nu + \\beta(y_{1j} - \\bar y_1)``.
Iteration 0 uses ``W = I`` and ``q = d`` (a linear start), then Eqn 15.17
is iterated with ``q`` and ``W`` of Eqns 15.18–15.19. `zero1`/`zero2`
constrain fixed-effect columns of the two traits.
"""
function joint_binary_quantitative(y1::AbstractVector, d::AbstractVector,
                                   X1::AbstractMatrix, Z1::AbstractMatrix,
                                   X2::AbstractMatrix, Z2::AbstractMatrix,
                                   Kinv::AbstractMatrix, G::AbstractMatrix, σ²e1::Real,
                                   r12::Real; ycentre = mean(y1), zero1 = Int[],
                                   zero2 = Int[], tol = 1e-8, maxiter = 50)
    n = length(y1)
    p1, q1, p2, q2 = size(X1, 2), size(Z1, 2), size(X2, 2), size(Z2, 2)
    β = r12 / sqrt(σ²e1) / sqrt(1 - r12^2)
    B = [1.0 0; -β 1]
    Gci = inv(B * G * B')
    K = Matrix{Float64}(Kinv)
    X1d, Z1d, X2d, Z2d = (Matrix{Float64}(M) for M in (X1, Z1, X2, Z2))
    ri = 1 / σ²e1
    ys = y1 .- ycentre
    τ, ν = zeros(p2), zeros(q2)
    o2 = p1 + q1                       # offset of the liability equations
    keep = setdiff(1:o2+p2+q2, [zero1; o2 .+ zero2])
    hist = Vector{Float64}[]
    b1, u1 = zeros(p1), zeros(q1)
    it = -1
    while it < maxiter
        it += 1
        if it == 0
            w, q = ones(n), Float64.(d)
        else
            μ = X2d * τ + Z2d * ν + β * ys
            q, w = zeros(n), zeros(n)
            for j in 1:n
                ϕ, F = _ϕ(μ[j]), _Φ(μ[j])
                d1, d2 = -ϕ / F, ϕ / (1 - F)
                q[j] = -(d[j] * d1 + (1 - d[j]) * d2)
                w[j] = μ[j] * q[j] + d[j] * d1^2 + (1 - d[j]) * d2^2
            end
        end
        WX, WZ = w .* X2d, w .* Z2d
        lhs = [ri*X1d'X1d ri*X1d'Z1d zeros(p1, p2) zeros(p1, q2)
               ri*Z1d'X1d ri*Z1d'Z1d+Gci[1, 1]*K zeros(q1, p2) Gci[1, 2]*K
               zeros(p2, p1) zeros(p2, q1) X2d'WX X2d'WZ
               zeros(q2, p1) Gci[2, 1]*K Z2d'WX Z2d'WZ+Gci[2, 2]*K]
        rhs = [ri * X1d'y1; ri * Z1d'y1 - Gci[1, 2] * (K * ν); X2d'q;
               Z2d'q - Gci[2, 2] * (K * ν)]
        x = zeros(length(rhs))
        x[keep] = Symmetric(lhs[keep, keep]) \ rhs[keep]
        b1, u1 = x[1:p1], x[p1+1:o2]
        Δτ, Δν = x[o2+1:o2+p2], x[o2+p2+1:end]
        τ .+= Δτ
        ν .+= Δν
        push!(hist, [b1; u1; τ; ν])
        it > 0 && maximum(abs, [Δτ; Δν]) < tol && break
    end
    JointBinaryResult(b1, u1, τ, ν, β, it, reduce(hcat, hist))
end
