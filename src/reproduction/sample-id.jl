"""
    function sampleID(mxy::BitMatrix, lmp::DataFrame, maf::Float64, ne::Int64)
Sample `ne` ID from the genotype `BitMatrix` `mxy`, and linkage map `lmp`. The
results are returned as `BitArray` `mxy` and `DataFrame` `lmp` with updated
frequencies. The loci are limited to `maf`.
"""
function sampleID(mxy::BitMatrix, lmp::DataFrame, maf::Float64, ne::Int64)
    0.0 ≤ maf ≤ 0.5 || error("maf must be between 0.0 and 0.5")
    nlc, nhp = size(mxy)
    nlc == nrow(lmp) || error("Genotype and map do not match")
    nid = nhp ÷ 2
    0 < ne ≤ nid || error("Improper number of ID to sample")

    id = sort(sample(1:nid, ne; replace = false))
    hps = repeat(id, inner = 2)
    hps .*= 2
    hps[1:2:end] .-= 1
    lmp.frq = vec(mean(mxy[:, hps], dims = 2))
    vlc = maf .< lmp.frq .< 1 - maf

    return copy(view(mxy, vlc, hps)), lmp[vlc, :]
end
