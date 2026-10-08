"""
    gibbs_sweep!(x, C, r; scale = 1.0, skip = Int[], rng = Random.default_rng()) -> x

One round of single-site Gibbs sampling of the location parameters of
mixed-model equations ``Cx = r`` (Mrode Eqns 18.11–18.12): in equation
order, ``x_j`` is drawn from

```math
N\\big((r_j - \\textstyle\\sum_{k \\ne j} c_{jk}x_k)/c_{jj},\\ \\text{scale}/c_{jj}\\big),
```

using the newest values of all other parameters. `scale` is
``\\sigma^2_e`` when `C` are the equations multiplied by the residual
variance (Mrode's form), `1` for unscaled equations. Equations in `skip`
are left unchanged (constraints). `C` must be symmetric; row `j` is read
from column `j`.
"""
function gibbs_sweep!(x::AbstractVector, C::SparseMatrixCSC, r::AbstractVector;
                      scale = 1.0, skip = Int[], rng::AbstractRNG = Random.default_rng())
    rv, nz, cp = rowvals(C), nonzeros(C), C.colptr
    sk = falses(length(x))
    sk[skip] .= true
    @inbounds for j in eachindex(x)
        sk[j] && continue
        s = r[j]
        d = 0.0
        for p in cp[j]:cp[j+1]-1
            i = rv[p]
            i == j ? (d = nz[p]) : (s -= nz[p] * x[i])
        end
        x[j] = s / d + randn(rng) * sqrt(scale / d)
    end
    x
end

"""
    GibbsResult

Output of [`gibbs`](@ref) and [`gibbs_mt`](@ref): posterior means of
the location parameters `x`, of the residual (co)variance `R` and of the
random-effect covariance matrices `G0`; the saved samples of the
variance parameters (one column per saved round; residual first, then
the lower triangles of each `G0`); and the number of saved samples.
"""
struct GibbsResult
    x::Vector{Float64}
    R::Matrix{Float64}
    G0::Vector{Matrix{Float64}}
    samples::Matrix{Float64}
    nsample::Int
end

# draw from an inverse Wishart / scaled inverse χ² with sum-of-squares matrix S
function _riw(rng, df, S)
    if size(S, 1) == 1
        return fill(S[1, 1] / rand(rng, Chisq(df)), 1, 1)
    end
    Matrix(rand(rng, InverseWishart(df, Matrix(Symmetric(S)))))
end

"""
    gibbs(y, X, effects; σ²e, νe = -2, Se = 0.0, νu = nothing, Vu = nothing,
          zero = Int[], niter = 10_000, burnin = 1_000, thin = 1,
          estimate_variances = true, rng = Random.default_rng()) -> GibbsResult

Gibbs sampler for variance components and breeding values of the
univariate mixed model of [`mme`](@ref) (Wang et al., 1993; Sorensen and
Gianola, 2002; Mrode Section 18.2): flat priors on fixed effects,
``u_k | G_{0k} \\sim N(0, G_{0k} \\otimes K_k)``, and

- ``\\sigma^2_e | \\cdot \\sim (e'e + \\nu_e S_e)/\\chi^2_{n + \\nu_e}``;
- ``G_{0k} | \\cdot \\sim IW(\\nu_{uk} + q_k,\\ U_k'K_k^{-1}U_k + V_{uk})``,
  a scaled inverse χ² for single-component effects.

The defaults (`νe = -2`, `Se = 0`, `νu = -(m+1)`, `Vu = 0`) are the
uniform priors of Mrode's examples. Locations are sampled with
[`gibbs_sweep!`](@ref) on the MME, which are rebuilt from the
precomputed ``W'W`` and ``W'y`` at every round. Samples after `burnin`
are kept every `thin` rounds. With `estimate_variances = false` the
variances stay at their starting values, and the posterior means of the
location parameters converge to the BLUP solutions.
"""
function gibbs(y::AbstractVector, X::AbstractMatrix, effects::AbstractVector{<:RandomEffect};
               σ²e::Real, νe::Real = -2, Se::Real = 0.0, νu = nothing, Vu = nothing,
               zero = Int[], niter::Integer = 10_000, burnin::Integer = 1_000,
               thin::Integer = 1, estimate_variances::Bool = true,
               rng::AbstractRNG = Random.default_rng())
    n = length(y)
    yv = Vector{Float64}(y)
    W = SparseMatrixCSC{Float64,Int}(hcat(sparse(X), (r.Z for r in effects)...))
    WW, Wy = W'W, W'yv
    p = size(X, 2)
    blocks = UnitRange{Int}[]
    o = p
    for r in effects
        push!(blocks, o+1:o+size(r.Z, 2))
        o += size(r.Z, 2)
    end
    ms = [ncomponents(r) for r in effects]
    νuv = νu === nothing ? [-(m + 1) for m in ms] : νu
    Vuv = Vu === nothing ? [zeros(m, m) for m in ms] : Vu
    Kis = [r.Kinv isa SparseMatrixCSC ? r.Kinv : sparse(r.Kinv) for r in effects]
    G0 = [copy(r.G0) for r in effects]
    ve = Float64(σ²e)
    x = zeros(size(W, 2))
    nsave = 0
    xs = zeros(length(x))
    Rs = 0.0
    Gs = [zeros(m, m) for m in ms]
    samples = Vector{Float64}[]
    for it in 1:niter
        pen = blockdiag(spzeros(p, p),
                        (kron(sparse(inv(Symmetric(G0[k]))), Kis[k]) for k in eachindex(effects))...)
        C = WW ./ ve + pen
        gibbs_sweep!(x, C, Wy ./ ve; skip = zero, rng = rng)
        if estimate_variances
            e = yv - W * x
            ve = (dot(e, e) + νe * Se) / rand(rng, Chisq(n + νe))
            for (k, r) in enumerate(effects)
                q = nlevels(r)
                U = reshape(view(x, blocks[k]), q, ms[k])
                S = U' * (Kis[k] * U) + Vuv[k]
                G0[k] = _riw(rng, νuv[k] + q, S)
            end
        end
        if it > burnin && (it - burnin) % thin == 0
            nsave += 1
            xs .+= x
            Rs += ve
            Gs .+= G0
            push!(samples, [ve; reduce(vcat, [[G[i, j] for (i, j) in _vech(size(G, 1))]
                                              for G in G0])])
        end
    end
    GibbsResult(xs ./ nsave, fill(Rs / nsave, 1, 1), Gs ./ nsave, reduce(hcat, samples),
                nsave)
end

"""
    gibbs_mt(Y, X, Z, Kinv; R0, G0, νe = -(t+1), Ve = 0, νu = -(t+1), Vu = 0,
             niter = 10_000, burnin = 1_000, thin = 1, estimate_variances = true,
             rng = Random.default_rng()) -> GibbsResult

Block Gibbs sampler for a multi-trait animal model with equal design
matrices (Jensen et al., 1994; Mrode Section 18.3):
``Y = XB + ZU + E``, rows of ``E`` ``\\sim N(0, R)``,
``\\operatorname{vec}(U) \\sim N(0, G \\otimes K)``. `Y` is `n × t` and may
contain `missing` or `NaN`; missing records are sampled from their
conditional distribution given the observed traits of the same record
(data augmentation) before each round.

Each round samples, for every fixed-effect level and every animal, the
`t` effects jointly (Eqns 18.23–18.24), then
``R \\sim IW(\\nu_e + n, E'E + V_e)`` and ``G \\sim IW(\\nu_u + q, U'K^{-1}U +
V_u)`` (Eqns 18.25–18.26). The defaults are uniform priors;
`estimate_variances = false` keeps `R0` and `G0` fixed. In the result,
`x` holds `vec([B; U])` (levels in rows, traits in columns).
"""
function gibbs_mt(Y::AbstractMatrix, X::AbstractMatrix, Z::AbstractMatrix,
                  Kinv::AbstractMatrix; R0::AbstractMatrix, G0::AbstractMatrix,
                  νe = nothing, Ve = nothing, νu = nothing, Vu = nothing,
                  niter::Integer = 10_000, burnin::Integer = 1_000, thin::Integer = 1,
                  estimate_variances::Bool = true, rng::AbstractRNG = Random.default_rng())
    n, t = size(Y)
    miss = map(v -> ismissing(v) || (v isa AbstractFloat && isnan(v)), Y)
    Yc = [miss[i, j] ? 0.0 : Float64(Y[i, j]) for i in 1:n, j in 1:t]
    Xs, Zs = sparse(Matrix{Float64}(X)), SparseMatrixCSC{Float64,Int}(Z)
    K = Kinv isa SparseMatrixCSC ? SparseMatrixCSC{Float64,Int}(Kinv) : sparse(Kinv)
    p, q = size(X, 2), size(Z, 2)
    νe_ = νe === nothing ? -(t + 1) : νe
    νu_ = νu === nothing ? -(t + 1) : νu
    Ve_ = Ve === nothing ? zeros(t, t) : Ve
    Vu_ = Vu === nothing ? zeros(t, t) : Vu
    R, G = Matrix{Float64}(R0), Matrix{Float64}(G0)
    B, U = zeros(p, t), zeros(q, t)
    E = Yc - Xs * B - Zs * U                         # residuals (updated in place)
    xx = vec(sum(abs2, Xs; dims = 1))
    zz = vec(sum(abs2, Zs; dims = 1))
    nsave = 0
    xs = zeros((p + q) * t)
    Rs, Gs = zeros(t, t), zeros(t, t)
    samples = Vector{Float64}[]
    for it in 1:niter
        # data augmentation of missing records
        for i in 1:n
            any(view(miss, i, :)) || continue
            o, m = findall(.!miss[i, :]), findall(miss[i, :])
            μ = Yc[i, :] - E[i, :]                    # current fitted values
            cμ = μ[m]
            cV = R[m, m]
            if !isempty(o)
                Kc = R[m, o] / R[o, o]
                cμ = cμ + Kc * (Yc[i, o] - μ[o])
                cV = cV - Kc * R[o, m]
            end
            new = cμ + cholesky(Symmetric(cV)).L * randn(rng, length(m))
            E[i, m] .+= new .- Yc[i, m]
            Yc[i, m] = new
        end
        Ri, Gi = inv(Symmetric(R)), inv(Symmetric(G))
        for k in 1:p                                   # fixed effects, Eqn 18.23
            rows, vals = rowvals(Xs)[nzrange(Xs, k)], nonzeros(Xs)[nzrange(Xs, k)]
            s = zeros(t)
            for (i, v) in zip(rows, vals)
                @views s .+= v .* (E[i, :] .+ v .* B[k, :])
            end
            Fd = cholesky(Symmetric(xx[k] .* Ri))      # D = U'U; U⁻¹z ~ N(0, D⁻¹)
            new = Fd \ (Ri * s) + Fd.U \ randn(rng, t)
            for (i, v) in zip(rows, vals)
                @views E[i, :] .-= v .* (new .- B[k, :])
            end
            B[k, :] = new
        end
        for j in 1:q                                   # animals, Eqn 18.24
            rows, vals = rowvals(Zs)[nzrange(Zs, j)], nonzeros(Zs)[nzrange(Zs, j)]
            s = zeros(t)
            for (i, v) in zip(rows, vals)
                @views s .+= v .* (E[i, :] .+ v .* U[j, :])
            end
            off = zeros(t)
            kjj = 0.0
            for p_ in nzrange(K, j)
                l = rowvals(K)[p_]
                l == j ? (kjj = nonzeros(K)[p_]) : (off .+= nonzeros(K)[p_] .* U[l, :])
            end
            Fd = cholesky(Symmetric(zz[j] .* Ri + kjj .* Gi))
            new = Fd \ (Ri * s - Gi * off) + Fd.U \ randn(rng, t)
            for (i, v) in zip(rows, vals)
                @views E[i, :] .-= v .* (new .- U[j, :])
            end
            U[j, :] = new
        end
        if estimate_variances
            R = _riw(rng, νe_ + n, E'E + Ve_)           # Eqn 18.25
            G = _riw(rng, νu_ + q, U' * (K * U) + Vu_)  # Eqn 18.26
        end
        if it > burnin && (it - burnin) % thin == 0
            nsave += 1
            xs .+= vec([B; U])
            Rs .+= R
            Gs .+= G
            push!(samples, [[R[i, j] for (i, j) in _vech(t)]; [G[i, j] for (i, j) in _vech(t)]])
        end
    end
    GibbsResult(xs ./ nsave, Rs ./ nsave, [Gs ./ nsave], reduce(hcat, samples), nsave)
end
