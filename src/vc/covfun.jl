"""
    covariance_function(G, t; k = size(G, 1), V = nothing, tmin = minimum(t),
                        tmax = maximum(t)) -> NamedTuple

Fit a covariance function (Kirkpatrick et al., 1990) of Legendre
polynomials of order `0:k-1` to a covariance matrix `G` observed at
ages/times `t`:

```math
G \\approx \\Phi C \\Phi'.
```

- Full order (`k = length(t)`): ``C = \\Phi^{-1} G \\Phi'^{-1}`` (Mrode
  Eqn 10.9).
- Reduced order (`k < length(t)`): generalized least squares on the
  ``t(t+1)/2`` distinct elements of `G`,
  ``\\check c = (X_s' V^{-1} X_s)^{-1} X_s' V^{-1} \\tilde g``, where `V`
  is the sampling covariance of the lower-triangle elements (column-major,
  ``G_{11}, G_{21}, …``); `V = I` (ordinary least squares) if omitted.

Returns `C` (`k × k`), `Φ`, the fitted `Ĝ = ΦCΦ'`, the coefficient
matrix `T` of the covariance function in powers of standardized age
(``\\Phi = M\\Lambda``, ``T = \\Lambda C \\Lambda'``; full order only,
else `nothing`), and Kirkpatrick's goodness-of-fit statistic `χ²` with
`df` degrees of freedom (`χ² = 0` for a full fit).
"""
function covariance_function(G::AbstractMatrix, t::AbstractVector; k::Integer = size(G, 1),
                             V = nothing, tmin = minimum(t), tmax = maximum(t))
    n = length(t)
    size(G) == (n, n) || throw(DimensionMismatch("G must be $n × $n"))
    1 ≤ k ≤ n || throw(ArgumentError("order k must be in 1:$n"))
    Φ = legendre(t, k; tmin, tmax)
    if k == n
        C = Φ \ G / Φ'
        Λ = _legendre_coefficients(k)
        x = @. -1 + 2 * (t - tmin) / (tmax - tmin)
        M = [xi^j for xi in x, j in 0:k-1]
        T = Λ * C * Λ'
        return (C = Symmetric(C) |> Matrix, Φ = Φ, Ĝ = Φ * C * Φ', T = T, M = M, χ² = 0.0,
                df = 0)
    end
    lo = [(i, j) for j in 1:n for i in j:n]            # lower triangle, column-major
    pr = [(i, j) for j in 1:k for i in j:k]            # distinct coefficients
    g = [G[i, j] for (i, j) in lo]
    Xs = [p == q ? Φ[i, p] * Φ[j, q] : Φ[i, p] * Φ[j, q] + Φ[i, q] * Φ[j, p]
          for (i, j) in lo, (p, q) in pr]
    W = V === nothing ? I : inv(Symmetric(Matrix{Float64}(V)))
    c = (Xs' * W * Xs) \ (Xs' * W * g)
    C = zeros(k, k)
    for (v, (p, q)) in zip(c, pr)
        C[p, q] = C[q, p] = v
    end
    r = g - Xs * c
    (C = C, Φ = Φ, Ĝ = Φ * C * Φ', T = nothing, M = nothing,
     χ² = V === nothing ? NaN : r' * W * r, df = length(lo) - length(pr))
end

"""
    _legendre_coefficients(k) -> Matrix

The `k × k` matrix ``\\Lambda'`` of normalized Legendre polynomial
coefficients, so that row `j+1` holds the coefficients of
``\\sqrt{(2j+1)/2}P_j(x)`` on ``1, x, …, x^{k-1}`` (Mrode Appendix G
uses ``\\Lambda`` with these as columns).
"""
function _legendre_coefficients(k::Integer)
    P = zeros(k, k)
    P[1, 1] = 1
    k > 1 && (P[2, 2] = 1)
    for j in 2:k-1                       # (j)P_j = (2j-1) x P_{j-1} - (j-1) P_{j-2}
        P[j+1, 2:end] .+= (2j - 1) / j .* P[j, 1:end-1]
        P[j+1, :] .-= (j - 1) / j .* P[j-1, :]
    end
    for j in 1:k
        P[j, :] .*= sqrt((2(j - 1) + 1) / 2)
    end
    P'
end
