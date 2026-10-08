"""
    mendelian_precision(sire, dam; F = nothing) -> Vector{Float64}

Inverse of the Mendelian-sampling variance factor of each animal,
``b_i = 1 / (1 - \\tfrac14 \\sum_{\\text{known } p}(1 + F_p))``, i.e.
``2``, ``4/3`` or ``1`` with both, one or no parents known and no
inbreeding. `sire`/`dam` codes outside `1:N` (e.g. `0` or group codes)
count as unknown. `F` are the inbreeding coefficients of animals `1:N`.
"""
function mendelian_precision(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer};
                             F = nothing)
    N = length(sire)
    Fv = F === nothing ? zeros(N) : F
    b = Vector{Float64}(undef, N)
    for i in 1:N
        v = 1.0
        1 ≤ sire[i] ≤ N && (v -= 0.25 * (1 + Fv[sire[i]]))
        1 ≤ dam[i] ≤ N && (v -= 0.25 * (1 + Fv[dam[i]]))
        b[i] = 1 / v
    end
    b
end

"""
    mendelian_precision(ped::DataFrame; inbreeding = true) -> Vector{Float64}

As above, from a pedigree; inbreeding by `RelationshipMatrices.nrm_diag`.
"""
mendelian_precision(ped::DataFrame; inbreeding::Bool = true) =
    mendelian_precision(ped.sire, ped.dam; F = inbreeding ? nrm_diag(ped) .- 1 : nothing)

"""
    ebv_partition(sire, dam, α, nrec, yd, a; F = nothing) -> NamedTuple

Partition the BLUP breeding values `a` of a univariate (animal or
repeatability) model into parent average (PA), yield deviation (YD)
and progeny contribution (PC), and compute progeny (daughter) yield
deviations (VanRaden and Wiggans, 1991; Mrode Eqns 4.8–4.13):

```math
\\hat a_i = n_1 \\mathrm{PA} + n_2 \\mathrm{YD} + n_3 \\mathrm{PC}
          = w_1 \\mathrm{PA} + w_2 \\mathrm{YD} + w_3 \\mathrm{DYD}.
```

# Arguments
- `sire`, `dam`: parents of animals `1:N` (`0` = unknown);
- `α`: ``\\sigma^2_e / \\sigma^2_a``;
- `nrec`: number of records of each animal;
- `yd`: mean yield deviation of each animal – its records adjusted for
  all fixed and non-genetic random effects (`0` without records);
- `a`: breeding-value solutions of animals `1:N`;
- `F`: optional inbreeding coefficients.

# Returns
A named tuple of vectors `PA, YD, PC, n1, n2, n3, DYD, w1, w2, w3`.
`PC` and `DYD` are `NaN` for animals without progeny.
"""
function ebv_partition(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}, α::Real,
                       nrec::AbstractVector, yd::AbstractVector, a::AbstractVector;
                       F = nothing)
    N = length(sire)
    b = mendelian_precision(sire, dam; F = F)
    par(p) = 1 ≤ p ≤ N ? a[p] : 0.0
    PA = [(par(sire[i]) + par(dam[i])) / 2 for i in 1:N]
    # progeny sums
    sp = zeros(N)    # Σ u_prog α / 2  (numerator of n3)
    spc = zeros(N)   # Σ u_prog α / 2 ⋅ (2â_o - â_mate)
    sd = zeros(N)    # Σ u_prog n2prog
    sdd = zeros(N)   # Σ u_prog n2prog (2YD_o - â_mate)
    for o in 1:N, (p, m) in ((sire[o], dam[o]), (dam[o], sire[o]))
        1 ≤ p ≤ N || continue
        u = b[o] / 2                    # u_prog: 1 or 2/3
        n2o = nrec[o] / (nrec[o] + α * b[o])
        sp[p] += 0.5α * u
        spc[p] += 0.5α * u * (2a[o] - par(m))
        sd[p] += u * n2o
        sdd[p] += u * n2o * (2yd[o] - par(m))
    end
    num1 = α .* b
    num2 = Float64.(nrec)
    den = num1 .+ num2 .+ sp
    PC = [sp[i] > 0 ? spc[i] / sp[i] : NaN for i in 1:N]
    DYD = [sd[i] > 0 ? sdd[i] / sd[i] : NaN for i in 1:N]
    num3w = 0.5α .* sd
    denw = num1 .+ num2 .+ num3w
    (PA = PA, YD = Vector{Float64}(yd), PC = PC, n1 = num1 ./ den, n2 = num2 ./ den,
     n3 = sp ./ den, DYD = DYD, w1 = num1 ./ denw, w2 = num2 ./ denw, w3 = num3w ./ denw)
end

"""
    ebv_partition_mt(sire, dam, G0, ZRZ, yd, a; F = nothing) -> NamedTuple

Multivariate version of [`ebv_partition`](@ref) (Mrode and Swanson,
2004; Mrode Eqns 6.5–6.13). For `t` traits:

```math
\\hat a_i = W_1 \\mathrm{PA} + W_2 \\mathrm{YD} + W_3 \\mathrm{PC}
          = M_1 \\mathrm{PA} + M_2 \\mathrm{YD} + M_3 \\mathrm{DYD},
```

with `t × t` weight matrices summing to `I`.

# Arguments
- `sire`, `dam`: parents of animals `1:N`;
- `G0`: `t × t` additive genetic covariance matrix;
- `ZRZ`: vector of the `t × t` blocks ``Z_i'R^{-1}Z_i`` of each animal
  (zero matrices for animals without records);
- `yd`: `N × t` yield deviations ``(Z_i'R^{-1}Z_i)^{-1}Z_i'R^{-1}(y_i -
  X_i\\hat b)`` (any value where `ZRZ` is zero);
- `a`: `N × t` breeding values.

# Returns
Vectors (over animals) of `PA, PC, DYD` (`t`-vectors) and of the weight
matrices `W1, W2, W3, M1, M2, M3`.
"""
function ebv_partition_mt(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer},
                          G0::AbstractMatrix, ZRZ::AbstractVector{<:AbstractMatrix},
                          yd::AbstractMatrix, a::AbstractMatrix; F = nothing)
    N, t = size(a)
    Gi = inv(Symmetric(Matrix{Float64}(G0)))
    b = mendelian_precision(sire, dam; F = F)
    row(i) = 1 ≤ i ≤ N ? a[i, :] : zeros(t)
    PA = [(row(sire[i]) + row(dam[i])) / 2 for i in 1:N]
    sp = zeros(N)                                  # Σ α_prog
    spc = [zeros(t) for _ in 1:N]                  # Σ α_prog (2â_o − â_mate)
    sW = [zeros(t, t) for _ in 1:N]                # Σ W2prog α_prog
    sWd = [zeros(t) for _ in 1:N]                  # Σ W2prog α_prog (2YD_o − â_mate)
    for o in 1:N, (p, m) in ((sire[o], dam[o]), (dam[o], sire[o]))
        1 ≤ p ≤ N || continue
        u = b[o] / 2
        W2o = (ZRZ[o] + b[o] * Gi) \ ZRZ[o]
        sp[p] += u
        spc[p] .+= u .* (2a[o, :] .- row(m))
        sW[p] .+= u .* W2o
        sWd[p] .+= u .* (W2o * (2yd[o, :] .- row(m)))
    end
    W1, W2, W3, M1, M2, M3 = (Vector{Matrix{Float64}}(undef, N) for _ in 1:6)
    PC = Vector{Vector{Float64}}(undef, N)
    DYD = Vector{Vector{Float64}}(undef, N)
    for i in 1:N
        D = ZRZ[i] + (b[i] + sp[i] / 2) * Gi
        W1[i] = D \ (b[i] * Gi)
        W2[i] = D \ ZRZ[i]
        W3[i] = D \ (0.5sp[i] * Gi)
        PC[i] = sp[i] > 0 ? spc[i] / sp[i] : fill(NaN, t)
        H = 0.5 * Gi * sW[i]
        D2 = ZRZ[i] + b[i] * Gi + H
        M1[i] = D2 \ (b[i] * Gi)
        M2[i] = D2 \ ZRZ[i]
        M3[i] = D2 \ H
        DYD[i] = sp[i] > 0 ? sW[i] \ sWd[i] : fill(NaN, t)
    end
    (PA = PA, PC = PC, DYD = DYD, W1 = W1, W2 = W2, W3 = W3, M1 = M1, M2 = M2, M3 = M3)
end
