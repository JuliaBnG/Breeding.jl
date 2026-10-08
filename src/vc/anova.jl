"""
    henderson3(y, X, Z; A = I) -> NamedTuple

Henderson's (1953) Method 3 for a sire model ``y = Xb + Zs + e``,
``\\operatorname{var}(s) = A\\sigma^2_s`` (Mrode Sections 17.2–17.4). With
``S = I - X(X'X)^-X'``, the sums of squares

- fixed effects ``F = y'X(X'X)^-X'y``,
- sires corrected for fixed effects ``S_s = y'SZ(Z'SZ)^-Z'Sy``, df = rank(Z'SZ),
- residual ``R = y'y - F - S_s``,

are equated to their expectations ``E(R) = df_R\\sigma^2_e`` and
``E(S_s) = df_S\\sigma^2_e + \\operatorname{tr}(Z'SZA)\\sigma^2_s``.
Returns the sums of squares, degrees of freedom, `ZSZ` (whose diagonal
holds the effective numbers of daughters) and the estimates `σ²e`,
`σ²s`.
"""
function henderson3(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix; A = I)
    n = length(y)
    Xd, Zd = Matrix{Float64}(X), Matrix{Float64}(Z)
    Px = Xd * pinv(Xd'Xd) * Xd'
    S = I - Px
    ZSZ = Zd' * S * Zd
    F = y' * Px * y
    Ss = (Zd'S * y)' * pinv(ZSZ) * (Zd'S * y)
    R = y'y - F - Ss
    dfX, dfS = rank(Xd), rank(ZSZ)
    dfR = n - dfX - dfS
    σ²e = R / dfR
    σ²s = (Ss - dfS * σ²e) / tr(ZSZ * A)
    (; F, S = Ss, R, dfX, dfS, dfR, ZSZ, σ²e, σ²s)
end

"""
    sire_contrasts(y, X, Z, A) -> NamedTuple

Canonical analysis of variance of a sire model (Mrode Section 17.4):
``Q`` with ``Q'(Z'SZ)Q = I`` and ``Q'(Z'SZAZ'SZ)Q = W`` (diagonal)
transforms ``Z'Sy`` into `df_S` independent contrasts ``c = Q'Z'Sy``
whose squares have expectations ``\\sigma^2_e + w_i\\sigma^2_s``. Returns
`Q`, the weights `w`, the contrast sums of squares `ss = c.^2`, and the
residual sum of squares `R` with `dfR` degrees of freedom; feed them to
[`fit_mean_squares`](@ref).
"""
function sire_contrasts(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix,
                        A::AbstractMatrix)
    h = henderson3(y, X, Z; A = A)
    B = Symmetric(h.ZSZ)
    E = eigen(B)
    keep = E.values .> 1e-10 * maximum(E.values)
    L = E.vectors[:, keep] ./ sqrt.(E.values[keep])'      # L'BL = I
    M = Symmetric(L' * B * A * B * L)
    F = eigen(M)
    o = sortperm(F.values; rev = true)
    Q = L * F.vectors[:, o]
    S = I - X * pinv(X'X) * X'
    c = Q' * (Z' * (S * y))
    (; Q, w = F.values[o], ss = c .^ 2, R = h.R, dfR = h.dfR)
end

"""
    fit_mean_squares(ms, df, coef; weighted = true, tol = 1e-10, maxiter = 1000)
        -> (θ, V)

Fit variance components ``\\theta`` to mean squares `ms` (each with `df`
degrees of freedom) whose expectations are `coef * θ` (one row per
mean square). With `weighted = false` ordinary least squares is used;
otherwise generalized least squares is iterated with
``\\operatorname{var}(ms_i) = 2E(ms_i)^2/df_i`` evaluated at the current
estimates (Mrode Section 17.5). Returns the estimates and their
sampling covariance matrix ``(C'V^{-1}C)^{-1}``.
"""
function fit_mean_squares(ms::AbstractVector, df::AbstractVector, coef::AbstractMatrix;
                          weighted = true, tol = 1e-10, maxiter = 1000)
    C = Matrix{Float64}(coef)
    θ = C \ ms
    weighted || return θ, inv(C'C) * (sum(abs2, ms - C * θ) / max(1, length(ms) - size(C, 2)))
    V = Diagonal(2 .* (C * θ) .^ 2 ./ df)
    for _ in 1:maxiter
        V = Diagonal(2 .* (C * θ) .^ 2 ./ df)
        θn = (C' * (V \ C)) \ (C' * (V \ ms))
        done = maximum(abs, θn - θ) < tol
        θ = θn
        done && break
    end
    θ, inv(C' * (V \ C))
end
