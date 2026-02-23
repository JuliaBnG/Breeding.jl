"""
    function sampleLoci(mxy::BitMatrix, lmp::DataFrame,sst::SNPSet...,)
Sample loci sets defined in `SNPSet`s `sst`. Each set in `sst` contains the name
of the set, the number of loci to sample, the minor allele frequency (MAF), and
a boolean showing if this set is exclusive to other loci sets. The results are
returned as `mxy` (`BitArray`) matrix and `lmp` (`DataFrame`) linkage map.
"""
function sampleLoci(mxy::BitMatrix, lmp::DataFrame, sst::SNPSet...)
    nlc, _ = size(mxy)
    nlc == nrow(lmp) || error("Genotype and map do not match")
    begin
        lc = 200 < nlc ? 200 : nlc
        frq = vec(mean(mxy[1:lc, :], dims = 2))
        f2q = lmp.frq[1:lc]
        all(frq .≈ f2q) || error("Genotype and map do not match")
    end
    nms = names(lmp)
    for nm in getfield.(sst, :name)
        nm ∈ nms && error("Duplicate name: ", nm)
        push!(nms, nm)
    end
    asnp = 1:nlc   # all available SNPs
    ssnp = Int64[] # selected SNPs

    for s in sst
        vld = s.maf .< lmp.frq .< 1 - s.maf
        pool = asnp[vld]
        s.exclusive && setdiff!(pool, ssnp)
        length(pool) < s.nlc && error("Not enough loci to sample from")
        lmp[!, s.name] .= false
        x = sort(sample(pool, s.nlc; replace = false))
        append!(ssnp, x)
        lmp[x, s.name] .= true
    end
    sort!(ssnp)
    mxy = copy(view(mxy, ssnp, :))
    return mxy, lmp[ssnp, :]
end
