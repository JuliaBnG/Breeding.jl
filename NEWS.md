# Breeding.jl release notes

## v0.2.2 (unreleased)

- `pev_candidates`: the prediction error (co)variance matrix of a set of
  candidates, from the mixed-model equations (`pev_candidates(m, k,
  candidates)`, solving one right-hand side per candidate) or from the
  selection-index form (`pev_candidates(K, ref, candidates, X, λ)`,
  which needs no inverse of `K`). The result is the `Ω` that `rocs`
  takes.
- `TM1997` and `OC_fixsex2` factor the relationship matrix of the active
  set by Cholesky instead of forming its inverse. Contributions agree
  with v0.2.1 to rounding (about 1e-11 on test problems), not bit for
  bit, and a relationship matrix that is not positive definite now
  throws `PosDefException`.
- Works with BnGStructs 0.2.6 and 0.3. With 0.3, `mtpblup` takes its
  column names from `BnGStructs.trait_columns`.
- Package metadata (`license`, `description`, URLs) in `Project.toml`;
  the tests are split by topic and run Aqua.jl.

## v0.2.1

- Support Julia 1.10 (LTS) through the latest release (1.13): `julia`
  compat lowered from `1.12` to `1.10`, stdlib compat relaxed to `1`.
  CI now tests 1.10, 1.12 and the latest release.

## v0.2.0

### Breaking

- Forward breeding simulation has moved to `BnGsim`. Removed from
  `Breeding`: `gene_drop`, `sum_map`, `sumMap`, `sampleID`, `sampleLoci`,
  `eqtl`, `tbv!`, `phenotype!` (and the internal `dummy_traits` and
  reproduction helpers). `StatsBase` is no longer a dependency.
- The unexported `TM2024` wrapper (which called an undefined `DOSop`) was
  removed.

### New

- `rocs` / `robust_ocs` and `ROCSResult`: robust optimum contribution
  selection of Fogg et al. (2026), covering TS, RTS, OCS and ROCS.
- `TM1997` and `OC_fixsex2` are now exported and documented; `TM1997`
  also accepts `(u, sex, A, K)` vectors. `OC_fixsex2` no longer prints
  debugging output.
- `zpz_hap!` accepts `BnGStructs.Haplotype` and `BnGStructs.Genotype`.
