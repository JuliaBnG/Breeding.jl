"""
    center_genotypes(M; p = column means / 2) -> (Z, p, k)

Centre an `n × m` matrix of allele dosages (`0, 1, 2`; individuals in
rows) as ``Z = M - 2\\mathbf{1}p'`` (VanRaden, 2008). Returns `Z`, the
allele frequencies `p` and ``k = 2\\sum_j p_j(1 - p_j)``, the scale
relating SNP and breeding-value variances, ``\\sigma^2_g =
\\sigma^2_a / k``.
"""
function center_genotypes(M::AbstractMatrix; p = vec(mean(M; dims = 1)) ./ 2)
    Z = Matrix{Float64}(M) .- 2 .* p'
    Z, Vector{Float64}(p), 2sum(p .* (1 .- p))
end

"""
    snp_from_gblup(Z, Ginv, a; k) -> Vector{Float64}

Back-solve SNP effects from GBLUP breeding values (Strandén and
Garrick, 2009): ``\\hat g = D Z' (ZDZ')^{-1} \\hat a = Z' G^{-1}\\hat
a / k`` with ``D = I / k`` and ``G = ZZ'/k``. `Ginv` must be the
inverse actually used in GBLUP (e.g. of ``G + 0.01I``).
"""
snp_from_gblup(Z::AbstractMatrix, Ginv::AbstractMatrix, a::AbstractVector; k::Real) =
    Z' * (Ginv * a) ./ k

"""
    gblup_index(G, ref, yc, λ; w = nothing) -> (dgv, rel)

Selection-index form of GBLUP (VanRaden, 2008; Mrode Eqns 11.10–11.12
and Section 11.9). With `yc` the records of the reference animals
`ref` corrected for fixed effects and ``\\lambda = \\sigma^2_e /
\\sigma^2_a``:

```math
\\hat a = G_{\\cdot r}\\,(G_{rr} + \\lambda R)^{-1} y_c,
\\qquad B = G_{\\cdot r}(G_{rr} + \\lambda R)^{-1} G_{r \\cdot},
```

for all animals in `G`, where ``R = \\operatorname{diag}(1/w)`` for
weights `w` (e.g. EDC; `I` if omitted). Theoretical reliabilities are
``b_{ii}/g_{ii}``.
"""
function gblup_index(G::AbstractMatrix, ref::AbstractVector{<:Integer}, yc::AbstractVector,
                     λ::Real; w = nothing)
    V = Symmetric(G[ref, ref] + λ * (w === nothing ? I : Diagonal(1 ./ w)))
    F = cholesky(V)
    Gr = G[:, ref]
    dgv = Gr * (F \ yc)
    B = Gr * (F \ Gr')
    dgv, diag(B) ./ diag(G)
end

"""
    base_allele_frequencies(M, Ainv, genotyped; λ = 0.01) -> (p, μ, d)
    base_allele_frequencies(M, ped::DataFrame, genotyped; λ = 0.01)

Base-population allele frequencies and gene content of ungenotyped
animals by the linear method of Gengler et al. (2007): every SNP is
analysed as a trait, ``q = \\mathbf{1}\\mu + Wd + e``, ``\\operatorname{var}(d)
= A\\sigma^2_d``, ``\\lambda = \\sigma^2_e/\\sigma^2_d``. `M` (`n_g × m`)
holds the dosages of the animals `genotyped` (indices into the `N × N`
pedigree inverse `Ainv`; with a pedigree `DataFrame` it is computed by
`RelationshipMatrices.ainv`).

The coefficient matrix is common to all SNPs, so it is factorized
once and all `m` right-hand sides are solved together. Returns the base
frequencies ``p = \\hat\\mu / 2``, the means, and the `N × m` BLUPs
`d` (predicted gene content is ``\\hat\\mu + \\hat d``).
"""
function base_allele_frequencies(M::AbstractMatrix, Ainv::AbstractMatrix,
                                 genotyped::AbstractVector{<:Integer}; λ = 0.01)
    N = size(Ainv, 1)
    ng = size(M, 1)
    length(genotyped) == ng || throw(DimensionMismatch("one row of M per genotyped animal"))
    W = incidence(genotyped, N)
    X = sparse(ones(ng, 1))
    lhs = [X'X X'W; W'X W'W+λ*SparseMatrixCSC{Float64,Int}(Ainv)]
    rhs = [X'M; W'M]
    S = cholesky(Symmetric(lhs)) \ Matrix{Float64}(rhs)
    μ = vec(S[1, :])
    μ ./ 2, μ, S[2:end, :]
end

base_allele_frequencies(M::AbstractMatrix, ped::DataFrame,
                        genotyped::AbstractVector{<:Integer}; λ = 0.01) =
    base_allele_frequencies(M, ainv(ped), genotyped; λ = λ)

"""
    pseudo_snps(haps::AbstractMatrix, blocks) -> (P, alleles)

Recode phased haplotypes into pseudo-SNPs (counts of each haplotype
allele within each block; e.g. Teissier et al., 2020). `haps` is
`2n × m` with rows `2i-1` and `2i` the two haplotypes of individual
`i`; `blocks` is a vector of column ranges. Haplotype alleles of a
block are numbered in order of first appearance. Returns the `n ×
Σ alleles` matrix of counts (`0`, `1`, `2`) and, per block, the allele
sequences.
"""
function pseudo_snps(haps::AbstractMatrix, blocks::AbstractVector{<:AbstractUnitRange})
    nh = size(haps, 1)
    iseven(nh) || throw(ArgumentError("haps must have two rows per individual"))
    n = nh ÷ 2
    cols = Matrix{Int8}[]
    alleles = Vector{Vector{Vector{eltype(haps)}}}()
    for blk in blocks
        seen = Dict{Vector{eltype(haps)},Int}()
        order = Vector{Vector{eltype(haps)}}()
        code = Vector{Int}(undef, nh)
        for h in 1:nh
            a = haps[h, blk]
            code[h] = get!(seen, a) do
                push!(order, a)
                length(order)
            end
        end
        C = zeros(Int8, n, length(order))
        for h in 1:nh
            C[(h+1)÷2, code[h]] += 1
        end
        push!(cols, C)
        push!(alleles, order)
    end
    reduce(hcat, cols), alleles
end
