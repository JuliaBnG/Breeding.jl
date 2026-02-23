"""
    tbv!(ped::DataFrame, xy::BitMatrix, lmp::DataFrame, trts)
Simulate true breeding values (TBV) for traits in `trts` using the genotype `xy`
and linkage map `lmp`. The TBV effects are calculated based on the QTL effects
in `lmp` and the genotype in `xy`. The resulting TBV values are added to the
`ped` DataFrame with column names prefixed by "tbv_". If all traits use the same
QTL set, the QTL loci are extracted once and used for all traits. Otherwise,
each trait's QTL loci are extracted individually.
"""
function tbv!(ped::DataFrame, xy::BitMatrix, lmp::DataFrame, trts)
    if length(unique(getfield.(trts, :QTL))) == 1
        qtl = lmp[:, trts[1].QTL] # QTL loci
        Q = qmat(xy, qtl) # QTL genotypes
        for t in trts
            eqtl = lmp[qtl, t.name] # QTL effects
            ped[!, "tbv_"*t.name] = zeros(size(ped, 1))
            copyto!(ped[!, "tbv_"*t.name], fastmv_int8_f64(Q, eqtl))
        end
        return ped
    else
        for t in trts
            qtl = lmp[:, t.QTL] # QTL loci
            Q = qmat(xy, qtl) # QTL genotypes
            eqtl = lmp[qtl, t.name] # QTL effects
            ped[!, "tbv_"*t.name] = fastmv_int8_f64(Q, eqtl)
        end
        return ped
    end
end
