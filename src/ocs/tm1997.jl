"""
    OCSContributions <: AbstractVector{Float64}

The contributions returned by [`TM1997`](@ref) and [`OC_fixsex2`](@ref): a
read-only vector (index it, broadcast over it, compare it with `≈` as a
`Vector{Float64}`) that also records whether the constraint was met.

# Fields
- `c::Vector{Float64}`: the contributions (the values of the vector).
- `feasible::Bool`: `true` if the contributions satisfy ``\\mathbf{c}'\\mathbf{A}\\mathbf{c}
  \\le K`` (to a relative tolerance of ``10^{-8}``); `false` if ``K`` could not be met, in
  which case `c` are the minimum-coancestry contributions that were found.
- `cAc::Float64`: ``\\mathbf{c}'\\mathbf{A}\\mathbf{c}`` of the contributions on the scale
  that `K` constrains (each sex contributing ½), with the fixed females included for
  `OC_fixsex2`.

Read the flag with [`feasible`](@ref). No warning is issued when ``K`` is not met; a caller
that iterates (e.g. one OCS step per generation) decides what to do.
"""
struct OCSContributions <: AbstractVector{Float64}
    c::Vector{Float64}
    feasible::Bool
    cAc::Float64
end

Base.size(r::OCSContributions) = size(r.c)
Base.getindex(r::OCSContributions, i::Int) = r.c[i]
Base.IndexStyle(::Type{OCSContributions}) = IndexLinear()

"""
    feasible(r::OCSContributions) -> Bool

Whether the contributions `r` satisfy the coancestry constraint ``\\mathbf{c}'\\mathbf{A}\\mathbf{c} \\le K``
they were computed for. See [`OCSContributions`](@ref).
"""
feasible(r::OCSContributions) = r.feasible

_meets(cAc, K) = cAc <= K + 1e-8 * max(abs(K), 1.0)

"""
    TM1997(ped::DataFrame, A::AbstractMatrix{<:Real}, K::Real) -> OCSContributions
    TM1997(u::AbstractVector{<:Real}, sex::AbstractVector, A::AbstractMatrix{<:Real}, K::Real) -> OCSContributions

Optimum genetic contribution selection with a predefined constraint on average
relationship `K` (Meuwissen, 1997, *J. Anim. Sci.* 75:934–940).

Maximizes expected genetic gain ``\\mathbf{c}'\\mathbf{u}`` subject to:
- ``\\mathbf{c}'\\mathbf{A}\\mathbf{c} \\le K``, where ``K = 2k_t`` is the maximum
  allowed average relationship among selected parents;
- ``\\mathbf{Q}'\\mathbf{c} = [\\tfrac{1}{2}, \\tfrac{1}{2}]'`` (equal total contribution
  from males and females);
- ``\\mathbf{c} \\ge \\mathbf{0}``.

The Lagrangian solution is iterated, dropping the candidates with negative
contributions, until all remaining contributions are positive. Candidates
once dropped are not re-admitted, so the result is not guaranteed to be the
optimum of the quadratic programme.

# Arguments
- `ped::DataFrame`: Candidate table with columns `tmi` (EBV or total merit index) and `sex` (`1`/`true` for males, `0`/`false` for females).
- `u::AbstractVector{<:Real}`: Vector of candidate EBVs or index values.
- `sex::AbstractVector`: Candidate sexes (`1`/`true` for males, `0`/`false` for females).
- `A::AbstractMatrix{<:Real}`: Symmetric positive definite relationship matrix of the candidates (`n × n`).
- `K::Real`: Upper constraint on ``\\mathbf{c}'\\mathbf{A}\\mathbf{c}``.

# Returns
An [`OCSContributions`](@ref): the contributions `c` (length `n`, summing to
`0.5` in each sex, all non-negative) as a vector, with `feasible(c)` and
`c.cAc`. If ``K`` is below the coancestry the iteration can reach, the
minimum-coancestry contributions are returned instead,
``\\mathbf{c} = \\tfrac12 \\mathbf{A}^{-1}\\mathbf{Q}(\\mathbf{Q}'\\mathbf{A}^{-1}\\mathbf{Q})^{-1}\\mathbf{1}``
over all candidates, again dropping negative contributions until none is
left; `feasible(c)` then says whether they meet ``K`` (`false` when ``K`` is
below the minimum coancestry, which is the usual case; `true` if only the
candidates dropped earlier made ``K`` unreachable, and these contributions
meet it without maximising gain). No warning is issued.
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
        F = cholesky(Symmetric(Matrix{Float64}(A[id, id])))
        AiQ = F \ Q
        Aiu = F \ u
        QAQi = inv(Symmetric(Q' * AiQ))                 # 2 × 2
        Qu = Q' * Aiu
        f1 = dot(u, Aiu) - dot(Qu, QAQi * Qu) # numerator
        f2 = 4K - sum(QAQi) # denominator
        if f2 ≤ 0                                       # K below the reachable coancestry
            rc = _min_coancestry(sex, A)
            cAc = dot(rc, A * rc)
            return OCSContributions(rc, _meets(cAc, K), cAc)
        end
        f1 = max(f1, 1e-10)
        λ₀ = sqrt(f1 / f2)
        λ = QAQi * (Qu .- λ₀)
        c = (Aiu - AiQ * λ) / (2λ₀)
        ix = findall(>(0), c) # indices of ID of next round
        if length(ix) == length(id)
            rc[id] = c
            cAc = dot(rc, A * rc)
            return OCSContributions(rc, _meets(cAc, K), cAc)
        end
        id = id[ix]
    end
end

# Minimum c'Ac with ½ per sex and c ≥ 0, over all candidates: the
# unconstrained minimum ½A⁻¹Q(Q'A⁻¹Q)⁻¹1, dropping the candidates with
# negative contributions until none is left (no re-admission).
function _min_coancestry(sex, A)
    nid = length(sex)
    id, rc = collect(1:nid), zeros(Float64, nid)
    while true
        Q = Float64[sex[id] .== 1 sex[id] .== 0]
        F = cholesky(Symmetric(Matrix{Float64}(A[id, id])))
        AiQ = F \ Q
        c = 0.5 .* (AiQ * (inv(Symmetric(Q' * AiQ)) * ones(2)))
        ix = findall(>(0), c)
        if length(ix) == length(id)
            rc[id] = c
            return rc
        end
        id = id[ix]
    end
end

"""
    OC_fixsex2(EBV::AbstractVector{<:Real}, A::AbstractMatrix{<:Real}, K::Real, c2::AbstractVector{<:Real}) -> OCSContributions

Optimize genetic contributions on the male side (`n1 = length(EBV)` candidates) given
fixed female contributions `c2` (`n2 = length(c2)` selected females), subject to
overall group relationship constraint `K` (Meuwissen, 1997).

# Arguments
- `EBV::AbstractVector{<:Real}`: EBVs of the `n1` male candidates.
- `A::AbstractMatrix{<:Real}`: Combined `(n1 + n2) × (n1 + n2)` relationship matrix where `A[1:n1, 1:n1]` is the male candidate block and `A[(n1+1):(n1+n2), (n1+1):(n1+n2)]` is the selected female block.
- `K::Real`: Target upper bound on average relationship, ``\\mathbf{c}'\\mathbf{A}\\mathbf{c}`` of all `n1 + n2` animals with each sex contributing ½.
- `c2::AbstractVector{<:Real}`: Fixed relative contributions of the `n2` females (internally scaled so ``\\sum c_2 = 0.5``).

# Returns
An [`OCSContributions`](@ref): the contributions of the `n1` male candidates
(scaled to sum to `1.0`), with `feasible` and `cAc` (``\\mathbf{c}'\\mathbf{A}\\mathbf{c}``
of males and females together, each sex ½). If ``K`` cannot be met, the male
contributions that minimise the coancestry given the females are returned,
``\\mathbf{c}_1 = \\mathbf{A}_{11}^{-1}(\\mu\\mathbf{1} - \\mathbf{A}_{12}\\mathbf{c}_2)``
with ``\\mu`` such that ``\\mathbf{1}'\\mathbf{c}_1 = \\tfrac12``, dropping negative
contributions, with `feasible == false`. No warning is issued.
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
        F = cholesky(Symmetric(Matrix{Float64}(A[ind, ind])))
        AI1 = F \ ones(length(ind))                  # row sums of A⁻¹
        AIu = F \ u
        sumAI = sum(AI1)
        AIA12 = F \ Matrix{Float64}(A[ind, (n1+1):(n1+n2)])
        A22_11 = A[(n1+1):(n1+n2), (n1+1):(n1+n2)] - A[(n1+1):(n1+n2), ind] * AIA12
        denominat =
            4 * (
                K - dot(c2_scaled, A22_11 * c2_scaled) -
                ((sum(AIA12 * c2_scaled) + 0.5)^2) / sumAI
            )
        numerat = dot(u, AIu) - sum(AIu)^2 / sumAI
        if denominat <= 0.0
            # minimum coancestry given the females: c₁ = A₁₁⁻¹(μ1 − A₁₂c₂), 1'c₁ = ½
            AIa = F \ A12c2[ind]
            μ = (0.5 + sum(AIa)) / sumAI
            c = μ .* AI1 .- AIa
        else
            numerat = max(numerat, 1e-10)
            lamb0 = sqrt(numerat / denominat)
            AIv = AIu - 2 * lamb0 .* (F \ A12c2[ind])
            lamb = (sum(AIv) - lamb0) / sumAI
            c = (AIv .- lamb .* AI1) ./ (2 * lamb0)
        end
        ind2 = findall(>=(0.0), c)
        if length(ind2) == length(ind)
            break
        end
        ind = ind[ind2]
    end

    cc = zeros(Float64, n)
    cc[ind] = c
    call = [cc; c2_scaled]
    cAc = dot(call, A * call)
    return OCSContributions(2 .* cc, _meets(cAc, K), cAc)
end
