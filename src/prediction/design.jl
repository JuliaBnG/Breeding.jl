"""
    incidence(x::AbstractVector; levels = sort(unique(x))) -> (SparseMatrixCSC, levels)

Incidence (design) matrix of a categorical factor `x`, with one column
per level in `levels`. Row `i` has a single `1` in the column of
`x[i]`. Records whose value is not in `levels` (e.g., a level dropped
to impose a constraint) get an all-zero row.

Returns the `length(x) × length(levels)` sparse matrix and the levels in
column order.
"""
function incidence(x::AbstractVector; levels = sort(unique(x)))
    col = Dict(l => j for (j, l) in enumerate(levels))
    I, J = Int[], Int[]
    sizehint!(I, length(x))
    sizehint!(J, length(x))
    for (i, v) in enumerate(x)
        j = get(col, v, 0)
        j == 0 && continue
        push!(I, i)
        push!(J, j)
    end
    sparse(I, J, ones(length(I)), length(x), length(levels)), levels
end

"""
    incidence(idx::AbstractVector{<:Integer}, n::Integer; w = 1.0) -> SparseMatrixCSC

Incidence matrix of `length(idx)` records on `n` columns, where record
`i` points to column `idx[i]` with weight `w` (a scalar or a vector).
Zero indices give empty rows. This is the usual `Z` relating records to
animals coded `1:n`.
"""
function incidence(idx::AbstractVector{<:Integer}, n::Integer; w = 1.0)
    I, J, V = Int[], Int[], Float64[]
    for (i, j) in enumerate(idx)
        j == 0 && continue
        1 ≤ j ≤ n || throw(BoundsError("column $j outside 1:$n"))
        push!(I, i)
        push!(J, j)
        push!(V, w isa Number ? w : w[i])
    end
    sparse(I, J, V, length(idx), n)
end

"""
    legendre(t::AbstractVector, k::Integer; tmin = minimum(t), tmax = maximum(t),
             normalized = true) -> Matrix{Float64}

Evaluate Legendre polynomials of order `0:k-1` at times `t`, after
standardizing `t` to ``[-1, 1]`` with ``x = -1 + 2(t - t_{min}) /
(t_{max} - t_{min})``. With `normalized = true` the polynomials are
scaled by ``\\sqrt{(2j + 1)/2}`` (Kirkpatrick et al., 1990), which is
the matrix ``\\Phi`` of Mrode's Appendix G.

Returns a `length(t) × k` matrix.
"""
function legendre(t::AbstractVector, k::Integer; tmin = minimum(t), tmax = maximum(t),
                  normalized = true)
    n = length(t)
    Φ = Matrix{Float64}(undef, n, k)
    @inbounds for i in 1:n
        x = -1 + 2 * (t[i] - tmin) / (tmax - tmin)
        p₀, p₁ = 1.0, x
        Φ[i, 1] = p₀
        k > 1 && (Φ[i, 2] = p₁)
        for j in 2:k-1 # Bonnet's recursion, P_j of order j
            p₀, p₁ = p₁, ((2j - 1) * x * p₁ - (j - 1) * p₀) / j
            Φ[i, j+1] = p₁
        end
    end
    if normalized
        for j in 1:k
            Φ[:, j] .*= sqrt((2(j - 1) + 1) / 2)
        end
    end
    Φ
end

"""
    mt_blocks(obs::AbstractMatrix{Bool}, Ms::AbstractVector{<:AbstractMatrix})
        -> Vector{SparseMatrixCSC{Float64,Int}}

Stack per-trait design matrices for a multivariate model with missing
records. `obs` is `n × t` (`true` = trait `j` observed on record row
`i`) and `Ms[j]` is the `n × pⱼ` design of trait `j` on all `n` rows.

Observations are stacked trait-major (all observed records of trait 1,
then trait 2, …), the order of [`mt_stack`](@ref). The `j`th returned
matrix is `nobs × pⱼ` with non-zeros only in trait `j`'s rows, so that
`hcat(blocks...)` is the block-diagonal fixed-effect design and the
blocks are the components of a multi-trait [`RandomEffect`](@ref).
"""
function mt_blocks(obs::AbstractMatrix{Bool}, Ms::AbstractVector{<:AbstractMatrix})
    n, t = size(obs)
    length(Ms) == t || throw(DimensionMismatch("need one design matrix per trait"))
    nobs = count(obs)
    blocks = Vector{SparseMatrixCSC{Float64,Int}}(undef, t)
    offset = 0
    for j in 1:t
        size(Ms[j], 1) == n || throw(DimensionMismatch("trait $j design has wrong rows"))
        rows = findall(view(obs, :, j))
        M = sparse(Ms[j][rows, :])
        I, J, V = findnz(M)
        blocks[j] = sparse(I .+ offset, J, V, nobs, size(M, 2))
        offset += length(rows)
    end
    blocks
end

"""
    mt_stack(Y::AbstractMatrix) -> (y::Vector{Float64}, obs::BitMatrix)

Stack an `n × t` phenotype matrix with `missing` (or `NaN`) entries
trait-major into a vector of observed values. Returns the vector and the
observation pattern used by [`mt_blocks`](@ref) and
[`mt_rinv`](@ref).
"""
function mt_stack(Y::AbstractMatrix)
    obs = BitMatrix(map(v -> !(ismissing(v) || (v isa AbstractFloat && isnan(v))), Y))
    y = Float64[Y[i, j] for j in axes(Y, 2) for i in axes(Y, 1) if obs[i, j]]
    y, obs
end

"""
    mt_rinv(R0::AbstractMatrix, obs::AbstractMatrix{Bool}) -> SparseMatrixCSC

Inverse residual (co)variance matrix of trait-major stacked records
(see [`mt_stack`](@ref)). Record row `i` with observed traits `S`
contributes ``(R_0[S,S])^{-1}``, so missing traits are handled exactly.
Inverses are cached per missing pattern.
"""
function mt_rinv(R0::AbstractMatrix, obs::AbstractMatrix{Bool})
    n, t = size(obs)
    size(R0) == (t, t) || throw(DimensionMismatch("R0 must be $t × $t"))
    pos = zeros(Int, n, t) # position of (i, j) in the stacked vector
    k = 0
    for j in 1:t, i in 1:n
        obs[i, j] && (pos[i, j] = (k += 1))
    end
    cache = Dict{BitVector,Matrix{Float64}}()
    I, J, V = Int[], Int[], Float64[]
    for i in 1:n
        S = BitVector(view(obs, i, :))
        any(S) || continue
        Ri = get!(() -> inv(Symmetric(Matrix{Float64}(R0[S, S]))), cache, S)
        ts = findall(S)
        for (b, tb) in enumerate(ts), (a, ta) in enumerate(ts)
            push!(I, pos[i, ta])
            push!(J, pos[i, tb])
            push!(V, Ri[a, b])
        end
    end
    sparse(I, J, V, k, k)
end

"""
    associate_incidence(group, animal, n; w = 1.0) -> SparseMatrixCSC

Incidence matrix ``Z_S`` of associative (social, indirect) genetic
effects (Muir, 2005; Bijma et al., 2007): record `i`, made by
`animal[i]` in `group[i]` (pen, cage, …), receives the associative
effect of every other member of its group, with weight `w` (e.g. a
dilution factor ``1/(n-1)``). `n` is the number of animal equations.
"""
function associate_incidence(group::AbstractVector, animal::AbstractVector{<:Integer},
                             n::Integer; w = 1.0)
    members = Dict{eltype(group),Vector{Int}}()
    for (g, a) in zip(group, animal)
        push!(get!(members, g, Int[]), a)
    end
    I, J = Int[], Int[]
    for (i, (g, a)) in enumerate(zip(group, animal))
        for b in members[g]
            b == a && continue
            push!(I, i)
            push!(J, b)
        end
    end
    sparse(I, J, fill(Float64(w), length(I)), length(animal), n)
end

"""
    linear_spline(t::AbstractVector, knots::AbstractVector) -> Matrix{Float64}

Covariables of a linear spline with the given (sorted) `knots` for
random-regression models (Misztal, 2006): a time `t` between knots
``T_i`` and ``T_{i+1}`` gets ``(T_{i+1} - t)/(T_{i+1} - T_i)`` on knot
`i` and ``(t - T_i)/(T_{i+1} - T_i)`` on knot `i+1`; outside the knots,
Mrode's rules ``t/T_1`` and ``T_n/t`` are applied to the end knot.
Returns a `length(t) × length(knots)` matrix with at most two non-zeros
per row.
"""
function linear_spline(t::AbstractVector, knots::AbstractVector)
    n = length(knots)
    issorted(knots) || throw(ArgumentError("knots must be sorted"))
    F = zeros(length(t), n)
    for (r, x) in enumerate(t)
        if x ≤ knots[1]
            F[r, 1] = x / knots[1]
        elseif x ≥ knots[n]
            F[r, n] = knots[n] / x
        else
            i = searchsortedlast(knots, x)
            w = (x - knots[i]) / (knots[i+1] - knots[i])
            F[r, i] = 1 - w
            F[r, i+1] = w
        end
    end
    F
end
