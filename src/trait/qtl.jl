"""
    eqtl(xy::BitMatrix, lmp::DataFrame, t::BnGStructs.AbstractTrait)
Simulate QTL effects for a single trait `t` using the genotype `BitMatrix`
`xy` and linkage map `lmp` `DataFrame`. The QTL effects are Normalized such that
the trait genetic variance has mean 0 and standard deviation `t.σₐ`, with a
default maximum deviation `ϵ = 1e-6`. The `lmp` is updated with the QTL effects.
"""
function eqtl(xy::BitMatrix, lmp::DataFrame, t::BnGStructs.AbstractTrait)
    if_match(xy, lmp) || error("Genotype and map do not match")
    string(t.QTL) ∈ names(lmp) || error("Missing $(t.QTL) column in lmp")
    qtl = lmp[!, t.QTL]
    Q = xy[qtl, 1:2:end] + xy[qtl, 2:2:end] # QTL genotypes
    lmp[!, t.name] .= 0.0 # add additive effect after sample_loci
    q = randn(sum(qtl)) # corrected to use length(qtl)
    norm_mv(Q, q; σ = t.σₐ) # => genotypes
    copyto!(view(lmp, qtl, t.name), q)
    nothing
end

"""
    eqtl(xy::BitMatrix, lmp::DataFrame, R::Matrix{Float64}, ts::Vector)
Simulate QTL effects for multiple traits `ts` using the genotype `BitMatrix`
`xy` and linkage map `lmp` `DataFrame`. The QTL effects are Normalized such that
the trait genetic variance has mean 0 and standard deviation `t.σₐ`, with a
default maximum deviation `ϵ = 1e-6`. The genetic correlation among `ts` is
defined by the matrix `R`. The `lmp` is updated with the QTL effects.
"""
function eqtl(xy::BitMatrix, lmp::DataFrame, R::Matrix{Float64}, ts::Vector)
    # Check parameters
    @info "  - Sampling QTL effects"
    if_match(xy, lmp) || error("Genotype and map do not match")
    issymmetric(R) && size(R, 1) == length(ts) || error("R and traits do not match")
    all(ts[1].QTL .== getfield.(ts, :QTL)) || error("Traits are using different QTL set")
    string(ts[1].QTL) ∈ names(lmp) || error("Missing $(ts[1].QTL) column in lmp")
    qtl = lmp[!, ts[1].QTL]
    nt, nqtl = length(ts), sum(qtl)
    Q = qmat(xy, qtl) # QTL genotypes
    U = cholesky(R).U
    qs = randn(nqtl, nt) * U
    print(' '^8)
    for (i, t) in enumerate(ts)
        lmp[!, t.name] .= 0.0
        print(' ', t.name)
        norm_mv(Q, view(qs, :, i); σ = t.σₐ)
        copyto!(view(lmp, qtl, t.name), view(qs, :, i))
    end
    println()
end
