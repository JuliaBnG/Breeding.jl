"""
    BayesResult

Posterior means from [`bayes_snp`](@ref): fixed effects `b`, SNP
effects `g`, SNP variances `σ²g` (per SNP for BayesA/B, a common value
for BayesC/Cπ, repeated), residual variance `σ²e`, `π` (the
probability that a SNP has no effect), the posterior inclusion
frequency of each SNP and the number of samples averaged.
"""
struct BayesResult
    b::Vector{Float64}
    g::Vector{Float64}
    σ²g::Vector{Float64}
    σ²e::Float64
    π::Float64
    inclusion::Vector{Float64}
    nsample::Int
end

"""
    bayes_snp(y, X, Z; method = :A, ν = 4.012, S = nothing, σ²g = nothing,
              σ²e = nothing, π = 0.3, g0 = 0.0, niter = 10_000, burnin = 3_000,
              nmh = 20, rng = Random.default_rng()) -> BayesResult

Gibbs samplers with residual updating (Legarra and Misztal, 2008) for
the SNP model ``y = Xb + Zg + e`` (Meuwissen et al., 2001; Habier et
al., 2011; Mrode Section 11.7):

- `:A` – SNP-specific variances, ``\\sigma^2_{g_i} | g_i \\sim
  (S + g_i^2)/\\chi^2_{\\nu+1}``;
- `:B` – as A, but ``\\sigma^2_{g_i} = 0`` with probability `π`;
  ``(\\sigma^2_{g_i}, g_i)`` are sampled jointly with `nmh`
  Metropolis–Hastings steps using the prior as proposal;
- `:C` – common variance ``\\sigma^2_g | g \\sim (S + g'g)/\\chi^2_{\\nu+k}``
  over the `k` fitted SNPs; a SNP enters with probability
  ``1 - \\pi`` through its indicator's full conditional;
- `:Cpi` – as C, with ``\\pi \\sim \\mathrm{Beta}(m - k + 1, k + 1)``.

Flat priors are used for `b` and for ``\\sigma^2_e``
(``\\sigma^2_e = e'e/\\chi^2_{n-2}``). `S` is the scale, by default
``\\tilde\\sigma^2_g(\\nu - 2)/\\nu`` from the starting SNP variance
`σ²g`, which defaults to ``\\operatorname{var}(y)/2 / (2\\sum p_j q_j)``
with ``p`` from `Z` (taken as uncentred dosages). Starting values:
``b = (X'X)^{-1}X'y``, `g0` for every SNP, `σ²e` = var(y) / 2.

Columns of `Z` are processed in place; the work per iteration is
``O(nm)`` and allocation-free.
"""
function bayes_snp(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix;
                   method::Symbol = :A, ν::Real = 4.012, S = nothing, σ²g = nothing,
                   σ²e = nothing, π::Real = 0.3, g0::Real = 0.0, niter::Integer = 10_000,
                   burnin::Integer = 3_000, nmh::Integer = 20,
                   rng::AbstractRNG = Random.default_rng())
    method in (:A, :B, :C, :Cpi) || throw(ArgumentError("method must be :A, :B, :C or :Cpi"))
    n, m = size(Z)
    nf = size(X, 2)
    Xd = Matrix{Float64}(X)
    Zd = Matrix{Float64}(Z)
    p = vec(mean(Zd; dims = 1)) ./ 2
    σ²g0 = σ²g === nothing ? var(y) / 2 / (2sum(p .* (1 .- p))) : Float64(σ²g)
    Sg = S === nothing ? σ²g0 * (ν - 2) / ν : Float64(S)
    b = Xd \ y
    g = fill(Float64(g0), m)
    vg = fill(σ²g0, m)                       # per-SNP or common variance
    δ = trues(m)
    ve = σ²e === nothing ? var(y) / 2 : Float64(σ²e)
    πc = Float64(π)
    e = y - Xd * b - Zd * g
    xx = vec(sum(abs2, Xd; dims = 1))
    zz = vec(sum(abs2, Zd; dims = 1))
    sb, sg, svg, sinc = zeros(nf), zeros(m), zeros(m), zeros(m)
    sve = sπ = 0.0
    ns = 0
    χe = Chisq(n - 2)
    for it in 1:niter
        ve = dot(e, e) / rand(rng, χe)                    # Eqn 11.21
        for j in 1:nf                                      # Eqn 11.20
            xj = view(Xd, :, j)
            rhs = dot(xj, e) + xx[j] * b[j]
            new = rhs / xx[j] + randn(rng) * sqrt(ve / xx[j])
            axpy!(b[j] - new, xj, e)
            b[j] = new
        end
        k = 0
        for i in 1:m
            zi = view(Zd, :, i)
            r = dot(zi, e) + zz[i] * g[i]                  # z'y*
            fit = true
            if method == :A
                vg[i] = (Sg + g[i]^2) / rand(rng, Chisq(ν + 1))   # Eqn 11.22
            elseif method == :B
                v = vg[i]
                ll(v) = (V = zz[i]^2 * v + zz[i] * ve; -0.5 * (log(V) + r^2 / V))
                l1 = ll(v)
                for _ in 1:nmh
                    vnew = rand(rng) < 1 - πc ? Sg / rand(rng, Chisq(ν)) : 0.0
                    l2 = ll(vnew)
                    if rand(rng) < exp(l2 - l1)
                        v, l1 = vnew, l2
                    end
                end
                vg[i] = v
                fit = v > 0
            else                                            # :C and :Cpi
                V1 = zz[i]^2 * vg[i] + zz[i] * ve
                V0 = zz[i] * ve
                l1 = -0.5 * (log(V1) + r^2 / V1) + log(1 - πc)
                l0 = -0.5 * (log(V0) + r^2 / V0) + log(πc)
                fit = rand(rng) < 1 / (1 + exp(l0 - l1))
            end
            new = if fit                                    # Eqn 11.23
                c = zz[i] + ve / vg[i]
                r / c + randn(rng) * sqrt(ve / c)
            else
                0.0
            end
            axpy!(g[i] - new, zi, e)
            g[i] = new
            δ[i] = fit
            k += fit
        end
        if method in (:C, :Cpi)                             # Eqn 11.28
            fill!(vg, (Sg + dot(g, g)) / rand(rng, Chisq(ν + k)))
            method == :Cpi && (πc = rand(rng, Beta(m - k + 1, k + 1)))
        end
        if it > burnin
            ns += 1
            sb .+= b
            sg .+= g
            svg .+= vg
            sinc .+= δ
            sve += ve
            sπ += πc
        end
    end
    BayesResult(sb ./ ns, sg ./ ns, svg ./ ns, sve / ns, sπ / ns, sinc ./ ns, ns)
end
