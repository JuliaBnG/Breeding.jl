"""
    mtpblup(ped, trts, giv, G; kw...)

Perform multi-trait PBLUP to estimate breeding values for traits `trts` in
pedigree `ped` using the inverse of relationship matrix `giv` and the genetic
(co)variance matrix `G`. The estimated breeding values (EBV) are updated in `ped`
in-place.

This is a simplified version of multi-trait PBLUP that assumes no other fixed
effects, no maternal effects, and no permanent environmental effects.

## Arguments
- `ped` : DataFrame, pedigree data frame containing columns "id", "sire", "dam",
  "sex", and phenotypes for traits in `trts`. The EBV are updated in this
  DataFrame.
- `trts` : Vector{Trait}, vector of Trait objects, the order should be the same
  as in the genetic (co)variance matrix `G`.
- `giv` : SparseMatrixCSC, the inverse of the numerator relationship matrix
  calculated from `ped`.
- `G` : Matrix, genetic (co)variance matrix for all traits.

- `kw...` : passed on to [`solve_mme`](@ref). The default is
  `method = :pcg` with [`pcg`](@ref)'s defaults (Jacobi preconditioner,
  relative residual `tol = 1e-12`), warm-started from the current EBV
  columns (`missing` as `0`).

## Returns
`nothing`; the EBV are written to the columns "ebv_"*`trait.name` of
`ped`.

## Notes
- The function will find common non-missing phenotypes for all the traits. Then
  set up the mixed model equations and solve for EBV using the
  preconditioned conjugate gradient method. The EBV are updated in `ped`
  in-place.
- The column names for phenotypes should be "ft_"*`trait.name`, and for EBV are
  "ebv_"*`trait.name`. For example, if trait carcass weight is named "wc", the
  phenotype column should be "ft_wc" and the EBV column should be "ebv_wc".
- Only the general means are fitted as fixed effects.
- Residuals are uncorrelated among traits, with ``\\sigma^2_{e,t} =
  G_{tt}(1 - h^2_t)/h^2_t``; the equations are
  ``\\mathbf{W'R^{-1}W} + \\mathbf{G}^{-1} \\otimes \\mathbf{A}^{-1}``
  (see [`mme`](@ref)).
"""
function mtpblup(ped, trts, giv, G; kw...)
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
    Z = incidence(findall(obs), nid)
    Zs = [kron(sparse(I, Tt, Tt)[:, t], Z) for t in 1:Tt] # trait-major records
    Y = Float64.(vec(Matrix(ped[obs, fts])))
    h² = getfield.(trts, :h²)
    σ²e = (1 .- h²) ./ h² .* diag(G) # no residual covariances
    m = mme(Y, X, [RandomEffect(Zs, giv, G)]; Rinv = repeat(1 ./ σ²e; inner = nob))
    x0 = [zeros(Tt); Float64.(coalesce.(vec(Matrix(select(ped, ebv))), 0.0))]
    sol = solve_mme(m; merge((; method = :pcg, x0), kw)...)
    for (i, e) in enumerate(ebv)
        fra = Tt + (i - 1) * nid + 1
        til = Tt + i * nid
        copyto!(ped[!, e], sol[fra:til])
    end
end
