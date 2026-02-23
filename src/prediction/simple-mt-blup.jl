"""" mtpblup(ped, trts, giv, G)

Perform multi-trait PBLUP to estimate breeding values for traits `trts` in
pedigree `ped` using the inverse of relationship matrix `giv` and the genetic
correlation matrix `G`. The estimated breeding values (EBV) are updated in `ped`
in-place.

This is a simplified version of multi-trait PBLUP that assumes no other fixed
effects, no maternal effects, and no permanent environmental effects.

## Arguments
- `ped` : DataFrame, pedigree data frame containing columns "id", "sire", "dam",
  "sex", and phenotypes for traits in `trts`. The EBV are updated in this
  DataFrame.
- `trts` : Vector{Trait}, vector of Trait objects, the order should be the same
  as in the genetic (co)variance matrix `G`.
- `giv` : SparseMatrixCSC, the inverse ofnumerical relationship matrix
  calculated from `ped`.
- `G` : Matrix, genetic (co)variance matrix for all traits.

## Returns
- `intercepts` : Vector, vector of trait intercepts.
- `Gsol` : Matrix, matrix of EBV for all traits (columns) and individuals
  (rows).

## Notes
- The function will find common non-missing phenotypes for all the traits. Then
  set up the mixed model equations and solve for EBV using Conjugate Gradient
  (CG) method. The EBV are updated in `ped` in-place.
- The column names for phenotypes should be "ft_"*`trait.name`, and for EBV are
  "ebv_"*`trait.name`. For example, if trait carcass weight is named "wc", the
  phenotype column should be "ft_wc" and the EBV column should be "ebv_wc".
- Only the general means are fitted as fixed effects.
"""
function mtpblup(ped, trts, giv, G)
    Tt = length(trts)  # total number of traits
    fts = "ft_" .* getfield.(trts, :name)
    ebv = "ebv_" .* getfield.(trts, :name)

    nid = size(ped, 1)
    obs = trues(nid)
    for f in fts
        obs .&= .!ismissing.(ped[!, f])
    end
    nob = sum(obs)
    X = sparse(kron(I(Tt), ones(nob)))
    Z₁ = begin
        Z = spzeros(nob, nid)
        for (i, j) in enumerate(findall(obs))
            Z[i, j] = 1
        end
        kron(I(Tt), Z)
    end
    Y = Float64.(vec(Matrix(ped[obs, fts])))
    h² = getfield.(trts, :h²)
    λ = (1 .- h²) ./ h²
    iG = inv(G) * Diagonal(λ)
    rhs = [
        X'Y
        Z₁'Y
    ]
    lhs = [ # assuming no permanent environmental effect
        X'X X'Z₁
        Z₁'X Z₁'Z₁ + kron(iG, giv)
    ]
    x0 = [zeros(Tt); vec(Matrix(select(ped, ebv)))]
    sol, _ = conjugate_gradient(lhs, rhs)
    for (i, e) in enumerate(ebv)
        fra = Tt + (i - 1) * nid + 1
        til = Tt + i * nid
        copyto!(ped[!, e], sol[fra:til])
    end
end
