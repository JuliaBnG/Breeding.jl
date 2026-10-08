"""
    AnimalModelIOD(y, fixed, animal, sire, dam, α; F = nothing)

Iteration-on-data (Schaeffer and Kennedy, 1986) representation of the
univariate animal model

```math
y = \\sum_f X_f b_f + Z a + e, \\quad \\operatorname{var}(a) = A\\sigma^2_a,
\\quad \\alpha = \\sigma^2_e / \\sigma^2_a,
```

that never forms the coefficient matrix. Only the data, the pedigree
and progeny lists are stored.

- `fixed`: vector of factors, each a vector of level indices `1:nₖ`
  (one per record);
- `animal`: animal index of each record;
- `sire`, `dam`: parents of animals `1:N`; `0` = unknown, `N + g` =
  unknown parent group `g` (phantom parents, QP-transformed
  ``A^{-1}``; Westell et al., 1988);
- `F`: optional inbreeding coefficients of animals `1:N`, used in the
  Mendelian-sampling precisions exactly as in
  `RelationshipMatrices.ainv`.

Equation order is: levels of each fixed factor, animals `1:N`, groups.
Use [`gauss_seidel`](@ref) on it for Gauss–Seidel iteration on data,
or [`pcg`](@ref) with the operator, `diag(m)` and `rhs(m)`.
"""
struct AnimalModelIOD
    y::Vector{Float64}
    fixed::Vector{Vector{Int}}
    nlev::Vector{Int}
    animal::Vector{Int}
    sire::Vector{Int}       # length nanim (= N + groups); groups have 0, 0
    dam::Vector{Int}
    δ::Vector{Float64}      # α / Mendelian sampling variance factor
    pptr::Vector{Int}       # progeny lists, CSR
    prog::Vector{Int}
    mate::Vector{Int}
    rptr::Vector{Int}       # records per animal, CSR
    rec::Vector{Int}
    foff::Vector{Int}       # equation offsets of the fixed factors
    α::Float64
end

function _csr(keys::AbstractVector{<:Integer}, n::Integer)
    cnt = zeros(Int, n + 1)
    for k in keys
        k > 0 && (cnt[k+1] += 1)
    end
    ptr = cumsum(cnt) .+ 1
    ptr[1] = 1
    ptr
end

function AnimalModelIOD(y::AbstractVector, fixed::AbstractVector{<:AbstractVector{<:Integer}},
                        animal::AbstractVector{<:Integer}, sire::AbstractVector{<:Integer},
                        dam::AbstractVector{<:Integer}, α::Real; F = nothing)
    N = length(sire)
    length(dam) == N || throw(DimensionMismatch("sire and dam differ in length"))
    nanim = max(N, maximum(sire), maximum(dam))
    Fv = F === nothing ? zeros(N) : Vector{Float64}(F)
    s = zeros(Int, nanim)
    d = zeros(Int, nanim)
    s[1:N] .= sire
    d[1:N] .= dam
    δ = zeros(nanim)
    for i in 1:N
        v = 1.0
        s[i] in 1:N && (v -= 0.25 * (1 + Fv[s[i]]))
        d[i] in 1:N && (v -= 0.25 * (1 + Fv[d[i]]))
        δ[i] = α / v
    end
    # progeny lists of every parent (animal or group)
    par = Int[]
    kid = Int[]
    mt = Int[]
    for i in 1:N
        s[i] > 0 && (push!(par, s[i]); push!(kid, i); push!(mt, d[i]))
        d[i] > 0 && (push!(par, d[i]); push!(kid, i); push!(mt, s[i]))
    end
    o = sortperm(par)
    pptr = _csr(par, nanim)
    rptr = _csr(animal, nanim)
    rec = sortperm(animal)
    nlev = [maximum(f) for f in fixed]
    foff = [0; cumsum(nlev)[1:end-1]]
    AnimalModelIOD(Vector{Float64}(y), [Vector{Int}(f) for f in fixed], nlev,
                   Vector{Int}(animal), s, d, δ, pptr, kid[o], mt[o], rptr, rec, foff,
                   Float64(α))
end

"""
    AnimalModelIOD(y, fixed, animal, ped::DataFrame, α; inbreeding = true)

Convenience constructor from a pedigree with columns `:sire` and `:dam`
(no unknown parent groups). With `inbreeding = true` the inbreeding
coefficients are computed by `RelationshipMatrices.nrm_diag`, so the
implicit ``A^{-1}`` equals `RelationshipMatrices.ainv(ped)`.
"""
function AnimalModelIOD(y::AbstractVector, fixed::AbstractVector{<:AbstractVector{<:Integer}},
                        animal::AbstractVector{<:Integer}, ped::DataFrame, α::Real;
                        inbreeding::Bool = true)
    F = inbreeding ? nrm_diag(ped) .- 1 : nothing
    AnimalModelIOD(y, fixed, animal, ped.sire, ped.dam, α; F = F)
end

neq(m::AnimalModelIOD) = sum(m.nlev) + length(m.δ)
_aoff(m::AnimalModelIOD) = sum(m.nlev)
_foff(m::AnimalModelIOD, f) = m.foff[f]

"""
    rhs(m::AnimalModelIOD) -> Vector{Float64}

Right-hand side ``[X'y; Z'y; 0]`` accumulated from the data.
"""
function rhs(m::AnimalModelIOD)
    r = zeros(neq(m))
    ao = _aoff(m)
    for (k, y) in enumerate(m.y)
        for f in eachindex(m.fixed)
            r[_foff(m, f)+m.fixed[f][k]] += y
        end
        r[ao+m.animal[k]] += y
    end
    r
end

"""
    diag(m::AnimalModelIOD) -> Vector{Float64}

Diagonal of the (implicit) coefficient matrix, used as the Jacobi
preconditioner of [`pcg`](@ref).
"""
function LinearAlgebra.diag(m::AnimalModelIOD)
    D = zeros(neq(m))
    ao = _aoff(m)
    for k in eachindex(m.y)
        for f in eachindex(m.fixed)
            D[_foff(m, f)+m.fixed[f][k]] += 1
        end
        D[ao+m.animal[k]] += 1
    end
    for i in eachindex(m.δ)
        δ = m.δ[i]
        δ == 0 && continue
        D[ao+i] += δ
        m.sire[i] > 0 && (D[ao+m.sire[i]] += δ / 4)
        m.dam[i] > 0 && (D[ao+m.dam[i]] += δ / 4)
    end
    D
end

"""
    (m::AnimalModelIOD)(v, d)

Overwrite `v` with ``\\mathbf{C}d``, reading the data and the pedigree
once – the matrix-free product of Section 19.6.
"""
function (m::AnimalModelIOD)(v::AbstractVector, d::AbstractVector)
    fill!(v, 0)
    ao = _aoff(m)
    @inbounds for k in eachindex(m.y)
        s = d[ao+m.animal[k]]
        for f in eachindex(m.fixed)
            s += d[_foff(m, f)+m.fixed[f][k]]
        end
        for f in eachindex(m.fixed)
            v[_foff(m, f)+m.fixed[f][k]] += s
        end
        v[ao+m.animal[k]] += s
    end
    @inbounds for i in eachindex(m.δ)
        δ = m.δ[i]
        δ == 0 && continue
        si, di = m.sire[i], m.dam[i]
        di_ = d[ao+i]
        ds = si > 0 ? d[ao+si] : 0.0
        dd = di > 0 ? d[ao+di] : 0.0
        v[ao+i] += δ * (di_ - 0.5ds - 0.5dd)
        si > 0 && (v[ao+si] += δ * (-0.5di_ + 0.25ds + 0.25dd))
        di > 0 && (v[ao+di] += δ * (-0.5di_ + 0.25ds + 0.25dd))
    end
    v
end

"""
    gauss_seidel(m::AnimalModelIOD; x0 = zeros, fix = Dict{Int,Float64}(),
                 tol = 1e-9, maxiter = 10_000, history = false)

Gauss–Seidel iteration on the data (Mrode Examples 19.3 and 19.4).
Fixed-effect levels are solved from records adjusted for the current
solutions of all other effects; each animal (and group) is solved from
its records, its parents and its progeny adjusted for their mates.
Equations in `fix` (equation index => value) are held at the given
values, e.g. to constrain group solutions.
"""
function gauss_seidel(m::AnimalModelIOD; x0 = zeros(neq(m)), fix = Dict{Int,Float64}(),
                      tol = 1e-9, maxiter = 10_000, history = false)
    x = Vector{Float64}(x0)
    for (k, v) in fix
        x[k] = v
    end
    xold = similar(x)
    ao = _aoff(m)
    nf = length(m.fixed)
    # records grouped by level for each fixed factor
    lvl = [(_csr(f, m.nlev[i]), sortperm(f)) for (i, f) in enumerate(m.fixed)]
    hist = history ? [copy(x)] : Vector{Float64}[]
    conv = Inf
    it = 0
    while it < maxiter
        it += 1
        copyto!(xold, x)
        for f in 1:nf
            ptr, rr = lvl[f]
            fo = _foff(m, f)
            for l in 1:m.nlev[f]
                haskey(fix, fo + l) && continue
                s = 0.0
                for p in ptr[l]:ptr[l+1]-1
                    k = rr[p]
                    adj = m.y[k] - x[ao+m.animal[k]]
                    for g in 1:nf
                        g == f && continue
                        adj -= x[_foff(m, g)+m.fixed[g][k]]
                    end
                    s += adj
                end
                n = ptr[l+1] - ptr[l]
                n > 0 && (x[fo+l] = s / n)
            end
        end
        for i in eachindex(m.δ)
            haskey(fix, ao + i) && continue
            D = 0.0
            R = 0.0
            for p in m.rptr[i]:m.rptr[i+1]-1         # own records
                k = m.rec[p]
                adj = m.y[k]
                for g in 1:nf
                    adj -= x[_foff(m, g)+m.fixed[g][k]]
                end
                R += adj
                D += 1
            end
            δ = m.δ[i]                                # parents
            if δ > 0
                D += δ
                m.sire[i] > 0 && (R += 0.5δ * x[ao+m.sire[i]])
                m.dam[i] > 0 && (R += 0.5δ * x[ao+m.dam[i]])
            end
            for p in m.pptr[i]:m.pptr[i+1]-1         # progeny
                o = m.prog[p]
                δo = m.δ[o]
                D += δo / 4
                mo = m.mate[p]
                R += 0.5δo * (x[ao+o] - (mo > 0 ? 0.5x[ao+mo] : 0.0))
            end
            D > 0 && (x[ao+i] = R / D)
        end
        history && push!(hist, copy(x))
        conv = _conv(x, xold)
        conv < tol && break
    end
    history ? (x, it, conv, reduce(hcat, hist)) : (x, it, conv)
end
