"""
    mtpblup(ped, trts, giv, G; kw...)

Perform multi-trait PBLUP to estimate breeding values for traits
`trts` in pedigree `ped` using the inverse relationship matrix `giv`
and genetic (co)variance matrix `G`. Estimated breeding values (EBV)
are written into `ped`.

This is a simplified multi-trait animal model:

```math
\\mathbf{y} = \\mathbf{Xb} + \\mathbf{Zu} + \\mathbf{e}, \\quad
\\operatorname{var}(\\mathbf{u}) = \\mathbf{G} \\otimes \\mathbf{A},
\\quad \\operatorname{var}(\\mathbf{e}) = \\mathbf{R} =
\\bigoplus_{t=1}^{T} \\sigma^2_{e,t} \\mathbf{I}
```

with trait-specific residual variances ``\\sigma^2_{e,t} = G_{tt}(1 -
h^2_t)/h^2_t`` and zero residual covariances. Only general trait means
are fitted as fixed effects; maternal and permanent environmental
effects are omitted.

## Arguments
- `ped` : `DataFrame` with phenotype columns `"ft_" * trait.name` for
  the traits in `trts`, one row per animal in the order of `giv`. EBV
  are written to `"ebv_" * trait.name`.
- `trts` : `Vector{Trait}`, ordered identically to rows and columns of
  `G`. Each trait must define `name` and `h²`.
- `giv` : `SparseMatrixCSC`, inverse additive relationship matrix
  ``\\mathbf{A}^{-1}`` of `ped`.
- `G` : `Matrix`, genetic (co)variance matrix among traits.
- `kw...` : passed to [`solve_mme`](@ref). Defaults to `method = :pcg`
  with [`pcg`](@ref)'s defaults (Jacobi preconditioner, relative
  residual `tol = 1e-12`), warm-started from the existing `"ebv_*"`
  column of each trait (`missing` treated as `0.0`); a trait without
  one starts from zeros.

## Returns
`nothing`; EBV are written to columns `"ebv_" * trait.name` of `ped`
(updating existing vectors in-place or creating new columns).

## Notes
- Only complete cases across all traits (records non-missing for every
  trait in `trts`) are included in the evaluation.
- The coefficient matrix is set up via [`mme`](@ref) as
  ``\\mathbf{W'R^{-1}W} + \\mathbf{G}^{-1} \\otimes \\mathbf{A}^{-1}``.
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
    nob > 0 ||
        throw(ArgumentError("no records with complete observations across all traits"))
    X = sparse(kron(I(Tt), ones(nob)))
    Z = incidence(findall(obs), nid)
    Zs = [kron(sparse(I, Tt, Tt)[:, t], Z) for t = 1:Tt] # trait-major records
    Y = Float64.(vec(Matrix(ped[obs, fts])))
    h² = getfield.(trts, :h²)
    σ²e = (1 .- h²) ./ h² .* diag(G) # no residual covariances
    m = mme(Y, X, [RandomEffect(Zs, giv, G)]; Rinv = repeat(1 ./ σ²e; inner = nob))
    x0 = zeros(Tt + Tt * nid)
    for (i, e) in enumerate(ebv) # warm start per trait, zeros if no column
        hasproperty(ped, e) || continue
        x0[(Tt+(i-1)*nid+1):(Tt+i*nid)] .= coalesce.(ped[!, e], 0.0)
    end
    sol = solve_mme(m; merge((; method = :pcg, x0), kw)...)
    for (i, e) in enumerate(ebv)
        fra = Tt + (i - 1) * nid + 1
        til = Tt + i * nid
        if hasproperty(ped, e)
            copyto!(ped[!, e], sol[fra:til])
        else
            ped[!, e] = sol[fra:til]
        end
    end
end
