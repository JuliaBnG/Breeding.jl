"""
    gamete(prt, hap, lms)

Generate a gamete `hap`, a **vector** view of a child's haplotype,
from `prt`, a view of a parents genotypes,
according crossovers generated from a linkage map summary, `lms`.
"""
function gamete(prt, hap, lms)
    cvs = crossovers(lms)
    h, l = cvs[1], 1            # starting
    for cv in cvs[2:end]
        copyto!(view(hap, l:cv), view(prt, l:cv, h))
        l = cv + 1
        h = 3 - h
    end
end

"""
    gamete(prt, hap, lms::DataFrame, tbp::Vector{UInt32})
Generate a gamete from a parent. `prt` contains the two haplotypes of a parent,
`hap` is the haplotype to be generated, `lms` is a DataFrame containing the
summary of the genome of the species. `tbp` is a vector of mutation positions in
the genome in base pairs. The function updates `hap`, and returns the genome
position of the crossover points in `gcs`. It can then be used to calculate the
true IBD segments. The crossover points for locus subscription are right
inclusive.
"""
function gamete(prt, hap, lms::DataFrame, tbp::Vector{UInt32})
    gcs, mcs = crossovers(lms, tbp)
    h, l = mcs[1], 1            # starting
    for cv in mcs[2:end]
        copyto!(view(hap, l:cv), view(prt, l:cv, h))
        l = cv + 1  # new segment starts after the crossover
        h = 3 - h   # switch haplotype
    end
    nlc = size(prt, 1)
    copyto!(view(hap, l:nlc), view(prt, l:nlc, h))  # copy remaining segment
    return gcs
end
