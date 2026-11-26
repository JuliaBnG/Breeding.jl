"""
    AnimalModel

A linear mixed model built from a table of records by
[`animal_model`](@ref). The equations are those of [`mme`](@ref) and
[`reml`](@ref):

- `y`: the observations, stacked trait-major (all records of trait 1,
  then trait 2, …), missing records left out;
- `X`: fixed-effect design, block-diagonal over traits, of full column
  rank;
- `effects`: the [`RandomEffect`](@ref)s, with their current `G0`;
- `residual`: [`Homogeneous`](@ref) for one trait, [`MultiTrait`](@ref)
  (one residual group per row of the table) for several;

and the maps back to the table:

- `traits`: trait names (columns of the table);
- `row`, `trait`: table row and trait index of every element of `y`;
- `fixed`: one row per fixed-effect level, with columns `trait`, `term`,
  `level` and `column` (column of `X`; `0` for a reference level, whose
  effect is set to zero);
- `levels`: the labels of the levels of every random effect, in the order
  of its equations;
- `components`: per random effect, the component names
  (`"<column>_<trait>"`), in the order of the rows of its `G0`.
"""
struct AnimalModel{R<:ResidualStructure}
    y::Vector{Float64}
    X::SparseMatrixCSC{Float64,Int}
    effects::Vector{RandomEffect}
    residual::R
    traits::Vector{Symbol}
    row::Vector{Int}
    trait::Vector{Int}
    fixed::DataFrame
    levels::Vector{Vector}
    components::Vector{Vector{String}}
end

Base.show(io::IO, m::AnimalModel) = print(io, "AnimalModel: ", length(m.traits), " trait(s), ",
    length(m.y), " records, ", size(m.X, 2), " fixed and ",
    sum(size(r.Z, 2) for r in m.effects; init = 0), " random equations")

"""
    animal_model(data; traits, fixed = Symbol[], covariates = Symbol[],
                 random = Pair[], ped = nothing, id = :id, sire = :sire, dam = :dam,
                 R0) -> AnimalModel

Build a linear mixed model from a table of records, without forming the
design matrices by hand:

```math
\\mathbf{y} = \\mathbf{Xb} + \\sum_k \\mathbf{Z}_k\\mathbf{u}_k + \\mathbf{e},
\\quad \\operatorname{var}(\\mathbf{u}_k) = \\mathbf{G}_{0k} \\otimes \\mathbf{K}_k,
\\quad \\operatorname{var}(\\mathbf{e}_i) = \\mathbf{R}_0[S_i, S_i],
```

where ``S_i`` are the traits recorded in row ``i`` of `data`. The result
goes to [`blup`](@ref), [`mme`](@ref) and [`reml`](@ref).

# Arguments
- `data`: any Tables.jl table, one row per record (repeated records of
  an animal are separate rows).
- `traits`: a column name, or a vector of them. A trait value that is
  `missing` or `NaN` is a missing record of that trait only.
- `fixed`: categorical fixed effects (column names), the same for every
  trait, or a `Dict` (or vector of pairs) `trait => columns`. Every trait
  has its own mean; the levels of a factor are those found among the
  records of the trait, sorted, and the first is the reference level,
  set to zero. Records of a trait with a `missing` fixed-effect value are
  dropped from that trait.
- `covariates`: numeric columns fitted as linear regressions within
  trait; given as `fixed`.
- `random`: a vector of `term => G0` or `term => (K, G0)`.
  - `term` is a column name, or a tuple of columns whose values index the
    same levels, e.g. `(:id, :dam)` for correlated direct and maternal
    genetic effects. A `missing` value gives no contribution for that
    record (e.g. an unknown dam).
  - `K` is `:I` (the default; levels independent, levels = the sorted
    values found in the table), `:A` (numerator relationships from
    `ped`; levels = the animals of the pedigree), or `(Kinv, labels)`, a
    `q × q` ``\\mathbf{K}^{-1}`` (e.g. ``\\mathbf{G}^{-1}``,
    ``\\mathbf{H}^{-1}``) whose rows belong to the animals `labels`.
  - `G0` is the covariance matrix among the components of the term,
    ordered column-major: `(column 1, trait 1), (column 1, trait 2), …,
    (column 2, trait 1), …`. With one trait and one column it may be a
    number.
- `ped`: the pedigree, a table with columns `id`, `sire` and `dam`
  holding labels of the same type as the table (unknown parents `0`,
  `""`, `"0"`, `missing` or `nothing`). It is recoded by
  `RelationshipMatrices.renumber`, which adds parents without a row as
  founders, and ``\\mathbf{A}^{-1}`` comes from `RelationshipMatrices.ainv`.
- `R0`: residual (co)variance matrix among traits (a number for one
  trait).

`id`, `sire` and `dam` name the pedigree columns. The starting values in
`G0` and `R0` are the covariances used by [`blup`](@ref) and the
starting point of [`reml`](@ref); `animal_model(m, r)` puts the REML
estimates `r` into the model.

Not covered: different random effects for different traits, unknown
parent groups, random regressions and heterogeneous residual classes,
which still need [`RandomEffect`](@ref) and [`mme`](@ref) directly.

# Examples
```julia
m = animal_model(df; traits = [:wwg, :pwg], fixed = [:herd, :sex],
                 random = [:id => (:A, G0)], ped = ped, R0 = R0)
r = reml(m)
est = blup(animal_model(m, r))
est.random[1]      # EBV of every animal of the pedigree, one column per trait
```
"""
function animal_model(data; traits, fixed = Symbol[], covariates = Symbol[],
                      random = Pair[], ped = nothing, id = :id, sire = :sire, dam = :dam, R0)
    df = DataFrame(data; copycols = false)
    trts = traits isa AbstractVector ? Symbol.(traits) : [Symbol(traits)]
    T = length(trts)
    T > 0 || throw(ArgumentError("give at least one trait"))
    allunique(trts) || throw(ArgumentError("a trait is listed twice"))
    fx = _per_trait(fixed, trts, "fixed")
    cv = _per_trait(covariates, trts, "covariates")
    for c in unique(vcat(trts, fx..., cv...))
        hasproperty(df, c) || throw(ArgumentError("the table has no column $c"))
    end
    R = R0 isa Real ? fill(Float64(R0), 1, 1) : Matrix{Float64}(R0)
    size(R) == (T, T) || throw(DimensionMismatch("R0 must be $T × $T"))

    # records: trait-major, rows of the table with the trait and its fixed effects known
    n = nrow(df)
    obs = falses(n, T)
    for (t, c) in enumerate(trts)
        col = df[!, c]
        need = [df[!, f] for f in vcat(fx[t], cv[t])]
        for i = 1:n
            v = col[i]
            _isknown(v) || continue
            v isa Real || throw(ArgumentError("trait $c has a non-numeric value $(repr(v))"))
            all(x -> !ismissing(x[i]), need) || continue
            obs[i, t] = true
        end
        any(view(obs, :, t)) || throw(ArgumentError("trait $c has no usable record"))
    end
    row = [i for t = 1:T for i = 1:n if obs[i, t]]
    trait = [t for t = 1:T for i = 1:n if obs[i, t]]
    ycols = [df[!, c] for c in trts]
    y = Float64[ycols[t][i] for (i, t) in zip(row, trait)]
    nobs = length(y)

    # fixed effects
    I_, J_, V_ = Int[], Int[], Float64[]
    ftab = DataFrame(trait = Symbol[], term = Symbol[], level = Any[], column = Int[])
    p = 0
    for t = 1:T
        recs = findall(==(t), trait)
        p += 1
        push!(ftab, (trts[t], :mean, missing, p))
        append!(I_, recs), append!(J_, fill(p, length(recs))), append!(V_, ones(length(recs)))
        for f in fx[t]
            vals = df[row[recs], f]
            lv = sort(unique(vals))
            push!(ftab, (trts[t], f, lv[1], 0))
            col = Dict(l => p + j - 1 for (j, l) in enumerate(lv) if j > 1)
            for (k, v) in zip(recs, vals)
                j = get(col, v, 0)
                j == 0 && continue
                push!(I_, k), push!(J_, j), push!(V_, 1.0)
            end
            for j = 2:length(lv)
                push!(ftab, (trts[t], f, lv[j], p + j - 1))
            end
            p += length(lv) - 1
        end
        for c in cv[t]
            p += 1
            push!(ftab, (trts[t], c, missing, p))
            vals = df[!, c]
            for k in recs
                v = vals[row[k]]
                v isa Real || throw(ArgumentError("covariate $c has a non-numeric value $(repr(v))"))
                push!(I_, k), push!(J_, p), push!(V_, Float64(v))
            end
        end
    end
    X = sparse(I_, J_, V_, nobs, p)
    rk = rank(qr(X))
    rk == p || throw(ArgumentError("the fixed effects are confounded: X has rank $rk " *
                                   "for $p columns (empty subclasses or collinear factors)"))

    # pedigree
    pedmap = nothing
    if ped !== nothing
        P = DataFrame(ped; copycols = false)
        for c in (id, sire, dam)
            hasproperty(P, c) || throw(ArgumentError("the pedigree has no column $c"))
        end
        r = RelationshipMatrices.renumber(P[!, id], P[!, sire], P[!, dam])
        pedmap = (labels = r.id, index = r.index, Ainv = SparseMatrixCSC{Float64,Int}(ainv(r.sire, r.dam)))
    end

    # random effects
    effects = RandomEffect[]
    levels = Vector[]
    comps = Vector{String}[]
    for spec in random
        spec isa Pair || throw(ArgumentError("a random term is `term => G0` or `term => (K, G0)`"))
        cols = spec.first isa Tuple ? Symbol.(collect(spec.first)) : [Symbol(spec.first)]
        K, G0 = spec.second isa Tuple ? spec.second : (:I, spec.second)
        for c in cols
            hasproperty(df, c) || throw(ArgumentError("the table has no column $c"))
        end
        m = length(cols) * T
        G = G0 isa Real ? fill(Float64(G0), 1, 1) : Matrix{Float64}(G0)
        size(G) == (m, m) || throw(DimensionMismatch(
            "G0 of $(Tuple(cols)) must be $m × $m ($(length(cols)) column(s) × $T trait(s))"))
        if K === :I
            lv = sort(unique(v for c in cols for v in df[row, c] if !ismissing(v)))
            index = Dict(l => j for (j, l) in enumerate(lv))
            Kinv = I
        elseif K === :A
            pedmap === nothing && throw(ArgumentError("K = :A needs a pedigree `ped`"))
            lv, index, Kinv = pedmap.labels, pedmap.index, pedmap.Ainv
        elseif K isa Tuple && length(K) == 2
            Kinv, lv = K
            size(Kinv) == (length(lv), length(lv)) ||
                throw(DimensionMismatch("Kinv must be $(length(lv)) × $(length(lv)), one row per label"))
            index = Dict(l => j for (j, l) in enumerate(lv))
            length(index) == length(lv) || throw(ArgumentError("duplicated labels of Kinv"))
        else
            throw(ArgumentError("K must be :I, :A or (Kinv, labels); got $(repr(K))"))
        end
        q = length(lv)
        Zs = SparseMatrixCSC{Float64,Int}[]
        names = String[]
        for c in cols
            vc = df[!, c]
            lvl = map(row) do i
                v = vc[i]
                ismissing(v) && return 0
                j = get(index, v, 0)
                j == 0 && throw(ArgumentError("$(repr(v)) in column $c is not a level of its " *
                                              (K === :A ? "pedigree" : "Kinv")))
                j
            end
            for t = 1:T
                push!(Zs, incidence([trait[k] == t ? lvl[k] : 0 for k = 1:nobs], q))
                push!(names, string(c, "_", trts[t]))
            end
        end
        push!(effects, RandomEffect(Zs, Kinv, G; name = join(cols, "+")))
        push!(levels, collect(lv))
        push!(comps, names)
    end

    res = T == 1 ? Homogeneous(R[1, 1]) : MultiTrait(trait, row, R)
    AnimalModel(y, X, effects, res, trts, row, trait, ftab, levels, comps)
end

_isknown(v) = !(ismissing(v) || v === nothing || (v isa AbstractFloat && isnan(v)))

function _per_trait(spec, trts, what)
    _cols(x) = x isa AbstractVector ? Symbol.(x) : [Symbol(x)]
    if spec isa AbstractDict || (spec isa AbstractVector && !isempty(spec) && all(x -> x isa Pair, spec))
        d = Dict(Symbol(k) => _cols(v) for (k, v) in spec)
        for k in keys(d)
            k in trts || throw(ArgumentError("$what given for $k, which is not a trait"))
        end
        [get(d, t, Symbol[]) for t in trts]
    else
        c = _cols(spec)
        [c for _ in trts]
    end
end

"""
    animal_model(m::AnimalModel, r::REMLResult) -> AnimalModel

The model `m` with the covariance matrices estimated by
[`reml`](@ref)`(m)`.
"""
function animal_model(m::AnimalModel, r::REMLResult)
    length(r.G0) == length(m.effects) || throw(ArgumentError("r does not come from this model"))
    effects = RandomEffect[RandomEffect(e.Z, e.Kinv, copy(G), e.name) for (e, G) in zip(m.effects, r.G0)]
    res = m.residual isa Homogeneous ? Homogeneous(r.σ²e) :
          MultiTrait(m.trait, m.row, r.R0[1])
    AnimalModel(m.y, m.X, effects, res, m.traits, m.row, m.trait, m.fixed, m.levels, m.components)
end

"""
    mme(m::AnimalModel) -> MME

The mixed-model equations of `m` at its current covariances: scaled by
``\\sigma^2_e`` for one trait, with the block-diagonal ``\\mathbf{R}^{-1}``
of the residual groups for several.
"""
function mme(m::AnimalModel)
    effs = [e for e in m.effects]
    m.residual isa Homogeneous && return mme(m.y, m.X, effs; σ²e = m.residual.σ²)
    mme(m.y, m.X, effs; Rinv = _rinv(m.residual))
end

# R⁻¹ of a MultiTrait structure with one class: block (R₀[S,S])⁻¹ per group
function _rinv(r::MultiTrait)
    length(r.R0) == 1 || throw(ArgumentError("heterogeneous residual classes are not supported here"))
    R0 = r.R0[1]
    groups = Dict{Int,Vector{Int}}()
    for (k, g) in enumerate(r.group)
        push!(get!(groups, g, Int[]), k)
    end
    cache = Dict{Vector{Int},Matrix{Float64}}()
    I_, J_, V_ = Int[], Int[], Float64[]
    for ks in values(groups)
        S = r.trait[ks]
        Ri = get!(() -> inv(Symmetric(R0[S, S])), cache, S)
        for (b, kb) in enumerate(ks), (a, ka) in enumerate(ks)
            push!(I_, ka), push!(J_, kb), push!(V_, Ri[a, b])
        end
    end
    n = length(r.trait)
    sparse(I_, J_, V_, n, n)
end

"""
    reml(m::AnimalModel; kw...) -> REMLResult

[`reml`](@ref) of the model `m`, started at its covariances; keywords
are passed on. Put the estimates into the model with
`animal_model(m, r)`.
"""
reml(m::AnimalModel; kw...) =
    reml(m.y, m.X, [e for e in m.effects]; residual = m.residual, kw...)

"""
    blup(m::AnimalModel; pev = false, kw...) -> (fixed, random, solutions)

Solve the mixed-model equations of `m` at its covariances (keywords go to
[`solve_mme`](@ref)) and return the estimates as tables:

- `fixed`: a `DataFrame` with columns `trait`, `term`, `level` and
  `estimate` (reference levels have `0.0`; the mean of a trait is the
  expectation of the reference subclass at covariates `0`);
- `random`: one `DataFrame` per random effect, a column `level` (the
  labels) and one column per component (`"<column>_<trait>"`); with
  `pev = true`, also `"pev_<column>_<trait>"`, the prediction-error
  variances from the selected inverse of the coefficient matrix;
- `solutions`: the solution vector of [`mme`](@ref)`(m)`.
"""
function blup(m::AnimalModel; pev::Bool = false, kw...)
    M = mme(m)
    x = solve_mme(M; kw...)
    b, us = solutions(M, x)
    fixed = select(m.fixed, :trait, :term, :level)
    fixed.estimate = [c == 0 ? 0.0 : b[c] for c in m.fixed.column]
    P = pev ? [Breeding.pev(M, k) for k in eachindex(m.effects)] : nothing
    random = map(eachindex(m.effects)) do k
        d = DataFrame(level = m.levels[k])
        for (j, c) in enumerate(m.components[k])
            d[!, c] = us[k][:, j]
        end
        if pev
            for (j, c) in enumerate(m.components[k])
                d[!, "pev_" * c] = P[k][:, j]
            end
        end
        d
    end
    (fixed = fixed, random = random, solutions = x)
end
