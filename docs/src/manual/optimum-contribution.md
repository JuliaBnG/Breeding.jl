# Optimum Contribution Selection

Optimum contribution selection (OCS) chooses the genetic contributions
$\mathbf{c}$ (or $\mathbf{w}$) of selection candidates to the next
generation so that expected genetic merit is maximised while the group
coancestry, and hence the rate of inbreeding, is controlled.
`Breeding.jl` provides two formulations:

- the **constraint** form of Meuwissen (1997): maximise $\mathbf{c}'\mathbf{u}$
  subject to an upper bound on $\mathbf{c}'\mathbf{A}\mathbf{c}$
  ([`TM1997`](@ref), [`OC_fixsex2`](@ref));
- the **penalty** form of Fogg et al. (2026), which also penalises the
  uncertainty of the EBVs ([`rocs`](@ref), alias `robust_ocs`).

The input EBVs, prediction-error (co)variances and relationship matrix
come directly from the rest of the package (e.g. [`solutions`](@ref),
[`pev_candidates`](@ref), [`pev`](@ref)) and from
`RelationshipMatrices.jl` (`nrm`, a genomic $\mathbf{G}$, …).

---

## Meuwissen (1997): constraint on average relationship

```math
\max_{\mathbf{c}} \ \mathbf{c}'\mathbf{u}
\quad \text{s.t.} \quad
\mathbf{c}'\mathbf{A}\mathbf{c} \le K, \quad
\mathbf{Q}'\mathbf{c} = \left[\tfrac12, \tfrac12\right]', \quad
\mathbf{c} \ge \mathbf{0}
```

where $\mathbf{Q}$ holds male and female indicators. With a target
inbreeding $k_t = k_{t-1} + \Delta F (1 - k_{t-1})$, $K = 2k_t$. The
Lagrangian solution is iterated, dropping candidates with negative
contributions, until all remaining contributions are positive. If $K$
is below the minimum attainable $\mathbf{c}'\mathbf{A}\mathbf{c}$, a
warning is issued and the minimum-coancestry contributions are
returned.

```julia
c = TM1997(DataFrame(tmi = ebv, sex = sex), A, K)   # sex: 1 male, 0 female
c = TM1997(ebv, sex, A, K)                          # equivalent vector method
```

[`OC_fixsex2`](@ref) solves the male side when the female contributions
are fixed (e.g. already-selected dams), with `A` ordered males first,
then the fixed females.

---

## Fogg et al. (2026): robust OCS

```math
\max_{\mathbf{l} \le \mathbf{w} \le \mathbf{u}}
\ \mathbf{w}'\bar{\boldsymbol{\mu}}
- \kappa \sqrt{\mathbf{w}'\boldsymbol{\Omega}\mathbf{w}}
- \frac{\lambda}{2} \mathbf{w}'\boldsymbol{\Sigma}\mathbf{w}
\quad \text{s.t.} \quad \mathbf{M}\mathbf{w} = \mathbf{m}
```

$\bar{\boldsymbol{\mu}}$ are the EBVs, $\boldsymbol{\Omega}$ their
prediction-error covariance matrix (or only its diagonal, the PEVs) and
$\boldsymbol{\Sigma}$ the relationship matrix. The two penalty weights
give four regimes:

| `λ` | `κ` | Regime |
|:---:|:---:|:-------|
| 0   | 0   | truncation selection (TS) |
| 0   | > 0 | robust truncation selection (RTS) |
| > 0 | 0   | optimum contribution selection (OCS) |
| > 0 | > 0 | robust OCS (ROCS) |

TS is solved exactly by sorting; OCS is a single convex quadratic
programme solved by an active-set method; RTS and ROCS iterate convex
quadratic subproblems whose fixed point satisfies the KKT conditions of
the robust problem. `sex` enforces $\sum w_i = \tfrac12$ within each sex
(otherwise $\sum w_i = 1$), and `lb`/`ub` bound individual contributions.

The scale of `λ` depends on whether $\boldsymbol{\Sigma}$ is a
relationship ($\mathbf{A}$, $\mathbf{G}$) or a coancestry
($\mathbf{A}/2$) matrix; `λ` and `κ` are trade-off weights, not
targets, and are usually swept to trace a Pareto frontier.

```julia
# Mrode Example 4.1: EBVs and their prediction-error covariance
m    = mme(y, X, [RandomEffect(Z, Ai, 20.0)]; σ²e = 40.0)
_, (u,) = solutions(m, solve_mme(m))
Ω    = pev_candidates(m, 1, 1:8)        # candidates: all eight animals
A    = nrm(ped)

# The same Ω from the selection-index form, for any relationship matrix
# K (A, G, H or IBD), records on the animals `ref`, and λ = σ²e/σ²a:
# Ω = pev_candidates(K, ref, candidates, X, 40.0 / 20.0; σ²a = 20.0)

r = rocs(vec(u), A; Ω = Ω, sex = sex, λ = 2.0, κ = 1.0, ub = 0.25)
r.w, r.merit, r.std_error, r.coancestry, r.converged
w, objective, merit, coancestry, std_error = r      # destructuring
```

---

## References

- Meuwissen, T. H. E. (1997). Maximizing the response of selection with
  a predefined rate of inbreeding. *J. Anim. Sci.* 75, 934–940.
- Fogg, Ortiz-Cuadros, Pocrnić, Hall and Gorjanc (2026). Robust optimal
  contribution selection. *Bioinformatics* 42(8), btag569.
