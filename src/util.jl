"""
    norm_mv(Q::AbstractMatrix, v; ϵ = 1e-6, μ=0., σ=1.)

Normalize QTL effect `v`, such that the `Q'v`, which are the trait genetic
values, have mean `μ` and standard deviation `σ`, with a default maximum
deviation `ε = 1e-6`.
"""
function norm_mv(Q::AbstractMatrix, v; ϵ = 1e-6, μ = 0.0, σ = 1.0)
    nqtl = size(Q, 1)
    tbv = Q'v
    m, s = mean(tbv), std(tbv)
    while abs(m - μ) > ϵ || abs(s - σ) > ϵ
        v .-= (m - μ) / nqtl
        v .*= σ / s
        tbv = Q'v
        m, s = mean(tbv), std(tbv)
    end
    tbv
end

"""
    if_match(xy::BitMatrix, lmp::DataFrame; slc = 200)
Check if the genotype `xy` and linkage map `lmp` match. The number of loci in
`xy` and `lmp` should be the same, and some random `slc = 200` loci of `xy`
should have the same allele frequencies as in `lmp`.
"""
function if_match(xy::BitMatrix, lmp::DataFrame; slc = 200)
    nlc = size(xy, 1)
    nlc == nrow(lmp) || return false
    slc < nlc && (slc = nlc)
    lc = sort(shuffle(1:nlc)[1:slc])
    frq = vec(mean(xy[lc, :], dims = 2))
    f2q = lmp.frq[lc]
    all(frq .≈ f2q) || return false
    true
end

"""
    qmat(xy::BitMatrix, qtl::Vector{Bool})
Calculate the QTL genotype matrix from the full genotype `xy` and the QTL
indicator `qtl`.
"""
function qmat(xy::BitMatrix, qtl::Vector{Bool})
    nid = size(xy, 2) ÷ 2
    Q = zeros(Int8, sum(qtl), nid)
    Threads.@threads for i = 1:nid
        Q[:, i] += xy[qtl, 2i-1]
        Q[:, i] += xy[qtl, 2i] # QTL genotypes
    end
    return Q
end

"""
This function calculates the transpose of 'a' multiplied by 'b' without creating
a large temporary matrix.

It is optimized with:
  - Multithreading (Threads.@threads) to run calculations in parallel.
  - @inbounds to remove array bounds checking.
  - @simd to encourage the compiler to use SIMD (Single Instruction, Multiple
    Data) instructions.
"""
function fastmv_int8_f64(mat::Matrix{Int8}, vec::Vector{Float64})
    m, n = size(mat)
    # Ensure dimensions are compatible for multiplication
    if length(vec) != m
        throw(DimensionMismatch("Matrix and vector dimensions are incompatible"))
    end

    # Pre-allocate the result vector. Use 'undef' because every element will be written to.
    result = Vector{Float64}(undef, n)

    # Parallelize the outer loop over the columns of the matrix.
    # Each thread will compute one or more dot products independently.
    Threads.@threads for i = 1:n
        # Use a local variable 'acc' to accumulate the sum for each dot product.
        # This is generally faster and safer for threading.
        acc = 0.0
        @inbounds @simd for j = 1:m
            # The core operation: multiply Int8 element by Float64 element.
            # The promotion to Float64 happens here, at the element level.
            acc += mat[j, i] * vec[j]
        end
        result[i] = acc
    end

    return result
end
