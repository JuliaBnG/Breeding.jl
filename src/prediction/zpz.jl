"""
    zpz!(ZZ, Z; accumulate = false, rows = 1024) -> ZZ

``\\mathbf{Z'Z}`` of an `n × m` genotype matrix `Z` (individuals in
rows, e.g. dosages `0, 1, 2` or centred genotypes), written into the
`m × m` matrix `ZZ`, or added to it with `accumulate = true`. Adding the
records of each new batch of individuals updates the SNP-BLUP
coefficient matrix without rebuilding it; keep ``\\lambda`` off the
diagonal of `ZZ` (add it to a copy) so that the result stays a pure
``\\mathbf{Z'Z}``.

For a BLAS element type of `ZZ` the lower triangle is computed with
`syrk`, converting `Z` in blocks of `rows` rows if its element type
differs (e.g. `Int8` dosages), and then mirrored to the upper
triangle. Integer dosages give exact sums, so accumulating over batches
gives the same `ZZ` as one call on all rows.
"""
function zpz!(
    ZZ::AbstractMatrix,
    Z::AbstractMatrix;
    accumulate::Bool = false,
    rows::Integer = 1024,
)
    rows ≥ 1 || throw(ArgumentError("rows must be positive, got $rows"))
    n, m = size(Z)
    size(ZZ) == (m, m) || throw(DimensionMismatch("ZZ must be $m × $m"))
    _zpz!(ZZ, Z, accumulate, rows)
    _mirror_lower!(ZZ)
end

function _zpz!(
    ZZ::StridedMatrix{T},
    Z::AbstractMatrix,
    accumulate,
    rows,
) where {T<:LinearAlgebra.BlasFloat}
    if Z isa StridedMatrix{T}
        BLAS.syrk!('L', 'T', one(T), Z, accumulate ? one(T) : zero(T), ZZ)
        return ZZ
    end
    n = size(Z, 1)
    accumulate || fill!(ZZ, zero(T))
    buf = Matrix{T}(undef, min(rows, n), size(Z, 2))
    for r = 1:rows:n
        k = min(rows, n - r + 1)
        B = view(buf, 1:k, :)
        copyto!(B, view(Z, r:(r+k-1), :))
        BLAS.syrk!('L', 'T', one(T), B, one(T), ZZ)
    end
    ZZ
end

function _zpz!(ZZ, Z, accumulate, rows) # e.g. integer ZZ; no overflow in eltype(Z)
    Zt = convert(AbstractMatrix{eltype(ZZ)}, Z)
    mul!(ZZ, Zt', Zt, true, accumulate)
end

function _mirror_lower!(ZZ)
    m = size(ZZ, 1)
    Threads.@threads for j = 2:m
        @inbounds for i = 1:(j-1)
            ZZ[i, j] = ZZ[j, i]
        end
    end
    ZZ
end

"""
    zpz_hap!(ZZ, H; accumulate = false, loci = :rows, block = 128) -> ZZ

``\\mathbf{Z'Z}`` of the genotype dosages ``z_{ki} = h_{k,i,1} +
h_{k,i,2}`` of phased 0/1 haplotypes `H`, computed by bit counting on
64 individuals per word without forming `Z`. `ZZ` (`m × m`) is
overwritten, or added to with `accumulate = true`; the counts are
integers, so accumulating over batches of individuals is exact (see
[`zpz!`](@ref)).

Layout of `H`:
- `loci = :rows` – `m × 2n`, columns `2k-1` and `2k` the haplotypes of
  individual `k` (the genotype layout of [`gene_drop`](@ref));
- `loci = :cols` – `n × 2m`, columns `2i-1` and `2i` the two haplotypes
  at locus `i`.

`H` may also be an integer matrix of 0/1 alleles; it is packed into a
`BitMatrix` first, and any other value (e.g. dosages `0, 1, 2`, which
belong in [`zpz!`](@ref)) throws an `ArgumentError`.

The lower triangle is computed in `block × block` tiles in parallel
over threads and mirrored to the upper one.
"""
function zpz_hap!(
    ZZ::AbstractMatrix,
    H::AbstractMatrix{Bool};
    accumulate::Bool = false,
    loci::Symbol = :rows,
    block::Integer = 128,
)
    block ≥ 1 || throw(ArgumentError("block must be positive, got $block"))
    loci in (:rows, :cols) || throw(ArgumentError("loci must be :rows or :cols"))
    r, c = size(H)
    iseven(c) || throw(
        ArgumentError(
            "H must have two columns per $(loci == :rows ? "individual" : "locus")",
        ),
    )
    m = loci == :rows ? r : c ÷ 2
    size(ZZ) == (m, m) || throw(DimensionMismatch("ZZ must be $m × $m"))

    # bit chunks of haplotype 1 and 2 at every locus, over individuals
    h1 = Vector{Vector{UInt64}}(undef, m)
    h2 = Vector{Vector{UInt64}}(undef, m)
    Threads.@threads for i = 1:m
        if loci == :rows
            h1[i] = BitVector(view(H, i, 1:2:c)).chunks
            h2[i] = BitVector(view(H, i, 2:2:c)).chunks
        else
            h1[i] = BitVector(view(H, :, 2i - 1)).chunks
            h2[i] = BitVector(view(H, :, 2i)).chunks
        end
    end

    nb = cld(m, block)
    tiles = [(ib, jb) for jb = 1:nb for ib = jb:nb]
    Threads.@threads :dynamic for t in eachindex(tiles)
        ib, jb = tiles[t]
        for j = ((jb-1)*block+1):min(jb*block, m)
            v1, v2 = h1[j], h2[j]
            for i = (ib==jb ? j : (ib-1)*block+1):min(ib*block, m)
                u1, u2 = h1[i], h2[i]
                s = 0
                @inbounds for w in eachindex(u1)
                    a, b = u1[w], u2[w]
                    x, y = v1[w], v2[w]
                    s +=
                        count_ones(a & x) +
                        count_ones(a & y) +
                        count_ones(b & x) +
                        count_ones(b & y)
                end
                @inbounds ZZ[i, j] = accumulate ? ZZ[i, j] + s : s
            end
        end
    end
    _mirror_lower!(ZZ)
end

function zpz_hap!(ZZ::AbstractMatrix, H::AbstractMatrix{<:Integer}; kw...)
    all(x -> x == 0 || x == 1, H) ||
        throw(ArgumentError("H must hold 0/1 alleles; use zpz! for dosages"))
    zpz_hap!(ZZ, BitMatrix(H .== 1); kw...)
end
