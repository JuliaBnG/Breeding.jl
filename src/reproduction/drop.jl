"""
    function drop(pg::AbstractArray, og::AbstractArray, pm, lms)
Drop haplotypes `pg` of parents into `og`, their offspring genotypes.
Parents of each offspring are defined in `pm`, which are rows of ``pa ma``.
Linkage map summary `lms` is from `summap`.

!!! ``Caution``:
- Merged data matrix from `MaCS` is `n-ID × n-loci`. Treat it with care.
- It is supposed all the genes to drop from are in `pg`.
- This function will be removed in the future.
"""
function drop(pg::AbstractArray, og::AbstractArray, pm, lms)
    nf = size(pm)[1]
    Threads.@threads for id = 1:nf
        ip = pm[id, 1]
        pa = view(pg, :, (2ip-1):2ip)
        zi = vec(view(og, :, 2id - 1))
        gamete(pa, zi, lms)
    end
    Threads.@threads for id = 1:nf
        im = pm[id, 2]
        ma = view(pg, :, (2im-1):2im)
        zi = vec(view(og, :, 2id))
        gamete(ma, zi, lms)
    end
end

"""
    gene_drop(
        pg::AbstractArray,
        og::AbstractArray,
        pm,
        lms::DataFrame,
        tbp::Vector{UInt32},
    )
Drop the recombined gametes in the parental generation `pg` into the offspring
generation `og` in place using the provided parental pairs `pm`. Crossovers are
generated using linkage map summary `lms` and mutation base pair position `tbp`
in the genome. The returned crossover point vectors have two leading values. The
first is reserved for batch, or generation, number. The second is the parent,
either `pa` or `ma`. These are the numbers in the previous batch, which are
provided in `pm`, or generation. They are then followed by the init haplotype
number, and the crossover points in the genome.
"""
function gene_drop(
    pg::AbstractArray,
    og::AbstractArray,
    pm,
    lms::DataFrame,
    tbp::Vector{UInt32},
)
    nf = size(pm)[1]
    cps = Vector{UInt32}[] # crossover points
    for (p, m) in eachrow(pm)
        push!(cps, [0, p])
        push!(cps, [0, m])
    end
    Threads.@threads for id = 1:nf
        ip, im = pm[id, :]
        # paternal
        pa = view(pg, :, (2ip-1):2ip)
        hp = vec(view(og, :, 2id - 1))
        gcsp = gamete(pa, hp, lms, tbp)
        append!(cps[2id-1], gcsp)
        # maternal
        ma = view(pg, :, (2im-1):2im)
        hm = vec(view(og, :, 2id))
        gcsm = gamete(ma, hm, lms, tbp)
        append!(cps[2id], gcsm)
    end
    df = DataFrame(
        generation = Int16[],
        parent = Int32[],
        ihap = Int8[],
        crossovers = Vector{UInt32}[],
    )
    for x in cps
        push!(df, (x[1], x[2], x[3], x[4:end]))
    end
    return df
end
