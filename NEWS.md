# Breeding.jl release notes

## v0.4.0 (2026-11-18)

### Breaking changes

- `TM1997` and `OC_fixsex2` return an `OCSContributions`, a read-only
  `AbstractVector{Float64}` of the contributions with two more fields:
  `feasible` (whether ``c'Ac \le K`` holds, read with `feasible(c)`) and
  `cAc` (the coancestry reached). Indexing, broadcasting (`2 .* c` is a
  `Vector{Float64}`) and `≈` work as before; code that needs a
  `Vector{Float64}` (e.g. to mutate it) uses `c.c` or `collect(c)`.
- No warning is issued when `K` cannot be met; the flag replaces it.

### Fixed

- `TM1997` with an unattainable `K` could return negative contributions:
  the minimum-coancestry fallback was the unconstrained solution on the
  candidates left after the last iteration. It is now computed over all
  candidates under ``c \ge 0`` (dropping negative contributions), and
  the tests check the KKT conditions of that minimum.
- `OC_fixsex2` with an unattainable `K` ignored the fixed females in its
  fallback (it returned equal row sums of ``A_{11}^{-1}``); it now
  returns the male contributions of minimum coancestry given the
  females, ``A_{11}^{-1}(\mu 1 - A_{12} c_2)``.

## v0.3.1 (2026-11-17)

- Compat: BnGStructs 0.4 (no change in Breeding; BnGStructs 0.4 removed
  `Species.nid` and `aTrait.QTL`, which Breeding does not use).

## v0.3.0 (2026-11-14)

### REML

- The REML rewrite of the `rr-reml` branch is merged:
  `reml(...; residual = ...)` with `Homogeneous`, `HetClasses`
  (heterogeneous residual classes) and `MultiTrait` residual structures,
  log-Cholesky AI steps, constraints (`fix`, `equal`), and a sparse
  back-end (`backend = :sparse`, the default above 2000 equations) that
  takes its traces from a selected inverse. `nested` builds nested
  design columns. New dependency: SelectedInversion.jl.
- `reml` no longer makes a dense copy of each `K⁻¹` to get `log|K|`. That
  copy alone needed 80 GB for 100 000 animals.

### Changed REML output

- `reml` records the log-likelihood of the final iterate as well, so
  `history` and `logLs` have the same number of columns (`logLs[k]`
  belongs to `history[:, k]`). In v0.2.1 `history` had one column more;
  code that paired them by dropping the last column of `history` must
  take `1:length(logLs)` instead. *(Added after release; it is also in
  the GitHub release notes.)*

### Prediction error variances

- `selected_inverse`: the elements of `C⁻¹` on the pattern of the Cholesky
  factor, at the cost of the factorisation.
- `pev` and `reliability` take the diagonal from the selected inverse
  instead of solving one column per level, and `reliability` no longer
  inverts `K⁻¹` densely. Values agree with the dense inverse (tests).

### Scale

`bench/scale.jl` times a single-trait animal model on a simulated
pedigree. On a Ryzen 9 3900X, 8 threads, Julia 1.13.3:

| Animals | all PEV | one REML iteration | PEV by one solve per animal (extrapolated from 20) |
| ---: | ---: | ---: | ---: |
| 10 000 | 0.19 s | 0.29 s | 18 s |
| 100 000 | 4.8 s | 9.2 s | about 3 h |

One run each, no replication.

Added after release: on a Ryzen 9 5950X with 93 GB, 16 threads, Julia
1.13.3 (one run, `bench/results/scale-2026-11-14-ryzen5950x.txt`):

| Animals | `ainv` | solve | all PEV | one REML iteration |
| ---: | ---: | ---: | ---: | ---: |
| 100 000 | 0.20 s | 1.5 s | 3.4 s | 6.6 s |
| 1 000 000 | 2.3 s | 71 s | 376 s | 892 s |

The process holding both runs peaked at 79 GB. From 10⁵ to 10⁶ animals
PEV took 110 times and REML 135 times as long, far more than linear; fill
in the Cholesky factor of this pedigree is the likely cause but has not
been measured. A first attempt at 10⁶ on a 31 GB desktop stalled and was
stopped, which is why `bench/scale.jl` now reports each stage on
`stderr`.

### Validation

`validation/wombat/validate.jl` compares `reml` with WOMBAT (version
26/05/2025) on simulated data: 3000 animals in 6 generations, 2500 with
records, 20 contemporary groups. For a univariate animal model and a
bivariate one with 19% of the records of the second trait missing, all
(co)variance estimates agree to a relative difference of 1.2e-5 or
less, the approximate standard errors to the four decimals WOMBAT
prints, and the maximised log-likelihoods to 4e-4. See
`validation/wombat/results.txt`.

## v0.2.2 (not registered; released as part of v0.3.0)

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
