"""
    phenotype!(ped::DataFrame, trts)
Generate phenotypes for traits `trts` in pedigree `ped`.
"""
function phenotype!(
    ped::DataFrame,
    trts, # one or more traits
)
    nid = nrow(ped)
    ft = Vector{Union{Missing,Float64}}(undef, nid) # continuous phenotypes
    tt = Vector{Union{Missing,Int8}}(undef, nid) # threshold phenotypes
    for t in trts
        σₑ = t.σₐ * sqrt((1 - t.h²) / t.h²)
        gt = ped[!, "tbv_"*t.name]
        for i = 1:nid
            ft[i] = gt[i] + rand(Normal(0, σₑ)) # phenotype without t.μ
        end
        if t.sex ≠ 2
            ft[ped.sex .== 1-t.sex] .= missing
        end
        if t isa BnGStructs.tTrait
            nt = length(t.threshold)
            for i = 1:nid
                ismissing(ft[i]) && continue
                n = findfirst(x -> x > ft[i], t.threshold)
                isnothing(n) ? (tt[i] = nt) : (tt[i] = n - 1)
            end
            ped[!, "ft_"*t.name] = tt
        else
            ped[!, "ft_"*t.name] = ft .+ t.μ
        end
    end
    ped
end
