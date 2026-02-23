"""
    sumMap(lmp::DataFrame; cM = 1e6)

Summary of a linkage map, which includes 3 columns, :chr, :pos(in bp), and :frq.
A DataFrame of 4-column for each chromosome is returned:
- numbering (1, 2, ...)
- length in Morgen
- number of loci
- beginning number of the first locus
"""
function sumMap(lmp::DataFrame; cM = 1e6)
    lms = DataFrame(
        chr = Int8[],
        len = Float64[], # as λ
        nlc = UInt32[],
        bgn = UInt32[],
    )
    bgn = 1
    for grp in groupby(lmp, :chr)
        chr = first(grp).chr
        len = (last(grp).pos - first(grp).pos) / (100cM)
        nlc = size(grp, 1)
        push!(lms, (chr, len, nlc, bgn))
        bgn += nlc
    end
    return lms
end

"""
    sum_map(chrLen::Vector{UInt32}; M = 1e8)
Summary of a linkage map, with information of chromosome length. The result
`DataFrame` has 3 columns:
1. `chr`: chromosome number
2. `len`, chromosome length in Morgan
3. `ebp`: ending base pair of the chromosome, in bp, of the genome
"""
function sum_map(chrLen::Vector{UInt32}; M = 1e8)
    nchr = length(chrLen)
    lms = DataFrame(chr = Int8.(1:nchr), len = chrLen / M, ebp = chrLen)
    for i = 2:nchr
        lms.ebp[i] += lms.ebp[i-1]
    end
    return lms
end
