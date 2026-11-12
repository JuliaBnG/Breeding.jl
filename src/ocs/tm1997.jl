"""
    TM1997(ped::DataFrame, A::AbstractMatrix{<:Real}, K::Real) -> Vector{Float64}
    TM1997(u::AbstractVector{<:Real}, sex::AbstractVector, A::AbstractMatrix{<:Real}, K::Real) -> Vector{Float64}

Optimum genetic contribution selection with a predefined constraint on average
relationship `K` (Meuwissen, 1997, *J. Anim. Sci.* 75:934–940).

Maximizes expected genetic gain ``\\mathbf{c}'\\mathbf{u}`` subject to:
- ``\\mathbf{c}'\\mathbf{A}\\mathbf{c} \\le K``, where ``K = 2k_t`` is the maximum
  allowed average relationship among selected parents;
- ``\\mathbf{Q}'\\mathbf{c} = [\\tfrac{1}{2}, \\tfrac{1}{2}]'`` (equal total contribution
  from males and females);
- ``\\mathbf{c} \\ge \\mathbf{0}``.

# Arguments
- `ped::DataFrame`: Candidate table with columns `tmi` (EBV or total merit index) and `sex` (`1`/`true` for males, `0`/`false` for females).
- `u::AbstractVector{<:Real}`: Vector of candidate EBVs or index values.
- `sex::AbstractVector`: Candidate sexes (`1`/`true` for males, `0`/`false` for females).
- `A::AbstractMatrix{<:Real}`: Symmetric positive definite relationship matrix of the candidates (`n × n`).
- `K::Real`: Upper constraint on ``\\mathbf{c}'\\mathbf{A}\\mathbf{c}``.

# Returns
- `Vector{Float64}`: Optimal parental contribution vector `c` of length `n` summing to `0.5` in each sex. If the relationship constraint `K` cannot be achieved, a warning is issued and the minimum-coancestry contributions are returned.
"""
function TM1997(ped::DataFrame, A::AbstractMatrix{<:Real}, K::Real)
    ["tmi", "sex"] ⊆ names(ped) || error("Not a proper pedigree")
    return TM1997(ped.tmi, ped.sex, A, K)
end

function TM1997(
    u_in::AbstractVector{<:Real},
    sex::AbstractVector,
    A::AbstractMatrix{<:Real},
    K::Real,
)
    issymmetric(A) || error("A is not symmetric")
    nid = length(u_in)
    length(sex) == nid || error("u and sex do not match")
    all(s -> s == 0 || s == 1, sex) ||
        throw(ArgumentError("sex must be coded 1/true (male) or 0/false (female)"))
    size(A, 1) == nid || error("A and ped do not match")
    id, rc = collect(1:nid), zeros(Float64, nid)

    while true
        u = Float64.(u_in[id])
        Q = Float64[sex[id] .== 1 sex[id] .== 0]
        Ai = inv(Matrix{Float64}(A[id, id]))
        QAQ = Q' * Ai * Q
        QAQi = inv(QAQ)
        f1 = u' * (Ai - Ai * Q * QAQi * Q' * Ai) * u # numerator
        f2 = 4K - sum(QAQi) # denominator
        if f2 ≤ 0
            rc[id] = 0.5 .* (Ai * Q * QAQi * ones(2)) # can't meet the constraint
            @warn "Cannot achieve constraint $K"
            return rc
        end
        f1 = max(f1, 1e-10)
        λ₀ = sqrt(f1 / f2)
        λ = QAQi * (Q' * Ai * u .- λ₀)
        c = Ai * (u - Q * λ) / (2λ₀)
        ix = findall(>(0), c) # indices of ID of next round
        if length(ix) == length(id)
            rc[id] = c
            return rc
        end
        id = id[ix]
    end
end

"""
    OC_fixsex2(EBV::AbstractVector{<:Real}, A::AbstractMatrix{<:Real}, K::Real, c2::AbstractVector{<:Real}) -> Vector{Float64}

Optimize genetic contributions on the male side (`n1 = length(EBV)` candidates) given
fixed female contributions `c2` (`n2 = length(c2)` selected females), subject to
overall group relationship constraint `K` (Meuwissen, 1997).

# Arguments
- `EBV::AbstractVector{<:Real}`: EBVs of the `n1` male candidates.
- `A::AbstractMatrix{<:Real}`: Combined `(n1 + n2) × (n1 + n2)` relationship matrix where `A[1:n1, 1:n1]` is the male candidate block and `A[(n1+1):(n1+n2), (n1+1):(n1+n2)]` is the selected female block.
- `K::Real`: Target upper bound on average relationship.
- `c2::AbstractVector{<:Real}`: Fixed relative contributions of the `n2` females (internally scaled so ``\\sum c_2 = 0.5``).

# Returns
- `Vector{Float64}`: Contributions of the `n1` male candidates (scaled to sum to `1.0`).
"""
function OC_fixsex2(
    EBV::AbstractVector{<:Real},
    A::AbstractMatrix{<:Real},
    K::Real,
    c2::AbstractVector{<:Real},
)
    c2_scaled = 0.5 .* Float64.(c2) ./ sum(c2)

    n = length(EBV)
    n1 = n
    n2 = length(c2_scaled)
    size(A) == (n1 + n2, n1 + n2) ||
        throw(DimensionMismatch("A must have size ($(n1 + n2), $(n1 + n2))"))
    ind = collect(1:n)
    A12c2 = A[1:n1, (n1+1):(n1+n2)] * c2_scaled
    c = zeros(Float64, n)
    for _ = 1:1000
        u = Float64.(EBV[ind])
        AI = inv(Matrix{Float64}(A[ind, ind]))
        sumAI = sum(AI)
        AIA12 = AI * A[ind, (n1+1):(n1+n2)]
        A22_11 = A[(n1+1):(n1+n2), (n1+1):(n1+n2)] - A[(n1+1):(n1+n2), ind] * AIA12
        denominat =
            4 * (
                K - dot(c2_scaled, A22_11 * c2_scaled) -
                ((sum(AIA12 * c2_scaled) + 0.5)^2) / sumAI
            )
        numerat = dot(u, AI * u) - sum(u' * AI)^2 / sumAI
        if denominat <= 0.0
            @warn "Cannot achieve constraint $K; minimizing relationships"
            c = vec(0.5 .* sum(AI; dims = 2) ./ sumAI)
        else
            numerat = max(numerat, 1e-10)
            lamb0 = sqrt(numerat / denominat)
            lamb = (sum(AI * (u - 2 * lamb0 .* A12c2[ind])) - lamb0) / sumAI
            c = AI * (u - 2 * lamb0 .* A12c2[ind] .- lamb) ./ (2 * lamb0)
        end
        ind2 = findall(>=(0.0), c)
        if length(ind2) == length(ind)
            break
        end
        ind = ind[ind2]
    end

    cc = zeros(Float64, n)
    cc[ind] = c
    return 2 .* cc
end
