"""
    crossovers(lms)

Give a linkage map summary `lms`, which can be from `sumMap`, this function
return a vector of crossover points along the whole genome.

DataFrame `lms` has four columns, chromosome number, its length in Morgen,
number of loci it contains, and the beginning number of its first locus in the
whole genome.

The first number is `1`, or `2`, the starting alternative haplotype. The vector
is then start from locus `1`, and loci that cross-over happens.  A
``crossover`` may also happen at the first locus of a chromsome. Otherwise, the
segment continues with the previous haplotype.  This eases the segment copy.  At
the beginning of a ``crossover`` happens on `rand([1:2])`.  The number of
cross-over occurred on a chromosome follows a Poisson distribution.

The crossover location is then sampled from a uniform distribution. It can can
also be U-shaped about the centromere, which might be implemented later.

The chromosome, or linkage group, should be in the order of the genotype file.
Or, unpredictable errors can happen.
"""
function crossovers(lms)
    pts = [Int32(rand(1:2))]           # crossover points
    for (_, λ, nlc, bgn) in eachrow(lms)
        bgn > 1 && rand(false:true) && push!(pts, bgn)
        nc = rand(Poisson(λ))
        append!(pts, sort(rand(1:nlc, nc) .+ (bgn - 1)))
    end
    push!(pts, last(lms).nlc - 1 + last(lms.bgn))
    pts
end

"""
    crossovers(lms::DataFrame, tbp::Vector{UInt32})
New version of `crossovers` function, which returns crossover bp in the genome, and
the mutation numbers. The former is for IBD calculation with Segment. The latter is
for the recombination to form offspring haplotype. Vector `tbp` is the pos column in
the linkage map `DataFrame`. `DataFrame` `lms` is just from the breed definition.
The first number is `1`, or `2`, the starting alternative haplotype in both vectors.
"""
function crossovers(lms::DataFrame, tbp::Vector{UInt32})
    gcvs = [UInt32(rand(1:2))]       # crossover points in genome bp
    mcps = [Int32(gcvs[1])]          # crossover with locus number
    bgn = UInt32(1)
    for (_, λ, ebp) in eachrow(lms)
        nc = rand(Poisson(λ))
        append!(gcvs, sort(rand(bgn:ebp, nc)))
        rand(false:true) && push!(gcvs, ebp)
        bgn = ebp + 1
    end
    gcvs[end] == lms.ebp[end] && pop!(gcvs)
    for bp in gcvs[2:end]
        i = searchsortedlast(tbp, bp)
        push!(mcps, i)
    end

    return gcvs, mcps
end
