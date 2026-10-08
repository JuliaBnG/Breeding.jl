"""
    ram_design(sire, dam, animal, σ²a, σ²e; F = nothing)
        -> (W, rinv, parents, m)

Design of the reduced animal model (Quaas and Pollak, 1980; Mrode
Section 4.5) for records on `animal` (indices into the pedigree
`sire`/`dam`, `0` = unknown parent).

Equations are set up for parents only. A parent's record maps to
itself; a non-parent's record maps to each known parent with ``1/2``
and carries the residual variance ``\\sigma^2_e + m_i\\sigma^2_a``, with
``m_i = 1 - \\tfrac14\\sum_{\\text{known } p}(1+F_p)`` the
Mendelian-sampling variance factor.

Returns the `nrec × np` sparse `W`, the inverse residual variances
`rinv` (use as `Rinv` in [`mme`](@ref) with a random effect on `W`
having ``K^{-1} = A_p^{-1}`` and variance `σ²a`), the sorted `parents`
(original indices; their sub-pedigree is closed, so `ainv` on it gives
``A_p^{-1}``), and the vector `m` of all animals.
"""
function ram_design(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer},
                    animal::AbstractVector{<:Integer}, σ²a::Real, σ²e::Real; F = nothing)
    N = length(sire)
    isparent = falses(N)
    for p in Iterators.flatten((sire, dam))
        p > 0 && (isparent[p] = true)
    end
    parents = findall(isparent)
    col = zeros(Int, N)
    col[parents] = 1:length(parents)
    m = 1 ./ mendelian_precision(sire, dam; F = F)
    I, J, V = Int[], Int[], Float64[]
    rinv = Vector{Float64}(undef, length(animal))
    for (k, i) in enumerate(animal)
        if isparent[i]
            push!(I, k); push!(J, col[i]); push!(V, 1.0)
            rinv[k] = 1 / σ²e
        else
            for p in (sire[i], dam[i])
                p > 0 && (push!(I, k); push!(J, col[p]); push!(V, 0.5))
            end
            rinv[k] = 1 / (σ²e + m[i] * σ²a)
        end
    end
    sparse(I, J, V, length(animal), length(parents)), rinv, parents, m
end

"""
    ram_backsolve(sire, dam, animal, yadj, ap, parents, σ²a, σ²e; F = nothing)
        -> Vector{Float64}

Breeding values of all animals from a reduced-animal-model solution:
parents take their solutions `ap` (ordered as `parents`); a non-parent
``i`` with ``n`` records gets ``\\mathrm{PA}_i + k\\,\\overline{(y -
x'\\hat b - \\mathrm{PA}_i)}`` with ``k = n / (n + \\sigma^2_e/(m_i
\\sigma^2_a))`` (Mrode Eqns 4.9 and 4.26). `yadj` are the records
already adjusted for the fixed effects. Animals with neither records
nor progeny get their parent average.
"""
function ram_backsolve(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer},
                       animal::AbstractVector{<:Integer}, yadj::AbstractVector,
                       ap::AbstractVector, parents::AbstractVector{<:Integer},
                       σ²a::Real, σ²e::Real; F = nothing)
    N = length(sire)
    a = fill(NaN, N)
    a[parents] = ap
    m = 1 ./ mendelian_precision(sire, dam; F = F)
    s = zeros(N)
    n = zeros(Int, N)
    for (k, i) in enumerate(animal)
        s[i] += yadj[k]
        n[i] += 1
    end
    for i in 1:N            # parents precede offspring
        isnan(a[i]) || continue
        pa = 0.5 * ((sire[i] > 0 ? a[sire[i]] : 0.0) + (dam[i] > 0 ? a[dam[i]] : 0.0))
        k = n[i] / (n[i] + σ²e / (m[i] * σ²a))
        a[i] = pa + (n[i] > 0 ? k * (s[i] / n[i] - pa) : 0.0)
    end
    a
end

"""
    ram_design(ped::DataFrame, animal, σ²a, σ²e; inbreeding = true)
        -> (W, rinv, parents, Apinv)

Reduced animal model from a pedigree: as the vector method, with
inbreeding from `RelationshipMatrices.nrm_diag` and, instead of the
Mendelian factors, the inverse relationship matrix of the parents
``A_p^{-1}`` = `RelationshipMatrices.ainv` of their (closed)
sub-pedigree, renumbered in the order of `parents`.
"""
function ram_design(ped::DataFrame, animal::AbstractVector{<:Integer}, σ²a::Real, σ²e::Real;
                    inbreeding::Bool = true)
    F = inbreeding ? nrm_diag(ped) .- 1 : nothing
    W, rinv, parents, _ = ram_design(ped.sire, ped.dam, animal, σ²a, σ²e; F = F)
    pos = zeros(Int, nrow(ped))
    pos[parents] = eachindex(parents)
    sub = DataFrame(sire = [s > 0 ? pos[s] : 0 for s in ped.sire[parents]],
                    dam = [d > 0 ? pos[d] : 0 for d in ped.dam[parents]])
    W, rinv, parents, ainv(sub)
end

"""
    ram_backsolve(ped::DataFrame, animal, yadj, ap, parents, σ²a, σ²e; inbreeding = true)

[`ram_backsolve`](@ref) with inbreeding from `RelationshipMatrices.nrm_diag`.
"""
ram_backsolve(ped::DataFrame, animal::AbstractVector{<:Integer}, yadj::AbstractVector,
              ap::AbstractVector, parents::AbstractVector{<:Integer}, σ²a::Real, σ²e::Real;
              inbreeding::Bool = true) =
    ram_backsolve(ped.sire, ped.dam, animal, yadj, ap, parents, σ²a, σ²e;
                  F = inbreeding ? nrm_diag(ped) .- 1 : nothing)
