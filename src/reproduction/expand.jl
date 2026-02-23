"""
    expand(
        xy::BitMatrix,
        lmp::DataFrame,
        lms::DataFrame,
        nₜ::Int64,
        t::Int64;
        maf = 0., # minor allele frequency threshold
    )
Expand the population defined in `xy` and `lmp` to `nₜ` individuals, in `t`
generations. The expansion is exponential, i.e., the number of individuals in
each `1:t` generation is `n₀⋅pⁱ``, with `n₀⋅pᵗ = nₜ`, where `n₀` is the number
of individuals in `xy`. The result genotype has a MAF, if `maf ∈ [0, 0.4)`. The
crossover points are returned in `cps`, which is a DataFrame with columns
`generation`, `sire`, `dam`, and `pos`. The crossovers in expansion can be used
to construct **A** matrix of the returned ID. The `lmp` DataFrame is updated
with the frequencies of the alleles in the returned genotypes.
"""
function expand(
    xy::BitMatrix,
    lmp::DataFrame,
    lms::DataFrame,
    nₜ::Int64,
    t::Int64;
    maf = 0.0, # minor allele frequency threshold
)
    nlc, nhp = size(xy)
    n₀ = nhp ÷ 2
    if_match(xy, lmp) || error("Genotype and map do not match")
    p = (nₜ / n₀)^(1 / t)
    prt, nₚ = copy(xy), n₀
    @info "  - Expanding $(nₚ) individuals to $(nₜ) in $(t) generations:"
    cps = DataFrame() # crossover points
    print(' '^8)
    for i = 1:t
        print(' ', i)
        nᵢ = round(Int, n₀ * p^i)
        off = BitArray(undef, nlc, 2nᵢ)
        pm = random_mate(rand(Bool, nₚ), nᵢ)
        cp = gene_drop(prt, off, pm, lms, lmp.pos)
        cp.generation .= Int16(i)
        append!(cps, cp)
        prt = nothing
        prt = off
        nₚ = nᵢ
    end
    print('\n')

    @info "  - Updating frequencies"
    nhₜ = size(prt, 2)
    frq = sum(prt, dims = 2) / nhₜ
    copyto!(lmp.frq, vec(frq))
    @info "  - Filtering by MAF threshold $maf"
    if 0.0 ≤ maf < 0.4
        vld = maf .< lmp.frq .< 1.0 - maf
    else
        vld = ones(Bool, nlc)
    end
    return prt[vld, :], lmp[vld, :], cps
end
