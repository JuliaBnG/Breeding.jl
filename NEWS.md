# Breeding.jl release notes

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
