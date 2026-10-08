# Breeding.jl

`Breeding.jl` is a Julia package for quantitative genetics: genetic
evaluation, genomic prediction, variance components and breeding-programme
simulation from pedigree and marker data.

[![Build Status](https://github.com/JuliaBnG/Breeding.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/Breeding.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-juliabng.github.io-blue.svg)](https://juliabng.github.io/#breeding)

The complete manual and API reference are available at
[juliabng.github.io/Breeding](https://juliabng.github.io/Breeding/).

---

## Features

- **Henderson's Mixed Model Equations (MME)**:
  Sparse coefficient matrix assembly for arbitrary fixed effects, multi-component random terms
  (direct + maternal, social/associative genetic effects, multi-trait, random regression with
  orthogonal Legendre polynomials or linear splines), prediction-error variances (PEV), and reliabilities ($r^2$).
- **Solvers & Iteration-on-Data (IOD)**:
  Sparse Cholesky, preconditioned conjugate gradient (`pcg`), Gauss–Seidel
  and SOR (`gauss_seidel` with `ω`), Jacobi (`jacobi`), and matrix-free
  iteration on data (`AnimalModelIOD`) supporting unknown parent groups.
- **Variance Component Estimation**:
  Restricted Maximum Likelihood (`reml`) via Average Information (AI) and Expectation-Maximization (EM),
  single-trait and multi-trait Gibbs samplers (`gibbs`, `gibbs_mt`), Henderson's Method 3 (`henderson3`),
  canonical sire contrasts, and covariance functions (`covariance_function`).
- **Genomic Prediction & Bayesian Alphabet**:
  VanRaden genotype centring (`center_genotypes`), GBLUP selection index (`gblup_index`), marker effect
  back-solving (`snp_from_gblup`), base allele frequency estimation via pedigree BLUP (`base_allele_frequencies`),
  phased pseudo-SNPs (`pseudo_snps`), and Bayesian alphabet regression (`bayes_snp`: BayesA, BayesB, BayesC, BayesCπ).
- **Threshold & Categorical Models**:
  Probit threshold models with Newton–Raphson Fisher scoring (`threshold_model`), category probability calculation
  (`category_probabilities`), and joint continuous-binary analysis (`joint_binary_quantitative`).
- **Pedigree Partitions & Reduced Models**:
  Mendelian sampling precisions (`mendelian_precision`), univariate and multivariate EBV partitioning
  (`ebv_partition`, `ebv_partition_mt`), and Reduced Animal Models (`ram_design`, `ram_backsolve`).
- **Forward Breeding Simulation**:
  Gene dropping with Poisson-distributed crossovers (`gene_drop`), founder and locus sampling
  (`sampleID`, `sampleLoci`), additive QTL effect simulation (`eqtl`), true breeding values (`tbv!`),
  and phenotype generation (`phenotype!`).

---

## Quick Start

Fit an animal model for weaning weight gain with fixed sex effect and random additive genetic effect:

```julia
using Breeding
using DataFrames
using RelationshipMatrices
using SparseArrays

# Pedigree and records
sire = [0, 0, 0, 1, 3, 1, 4, 3]
dam  = [0, 0, 0, 0, 2, 2, 5, 6]
calf = [4, 5, 6, 7, 8]
sex  = ["M", "F", "F", "M", "M"]
wwg  = [4.5, 2.9, 3.9, 3.5, 5.0]

# Setup design matrices
X, lv = incidence(sex; levels = ["M", "F"])
Z = incidence(calf, 8)

# Additive relationship matrix inverse A⁻¹
ped = DataFrame(sire = sire, dam = dam)
Ai = ainv(ped)

# Setup random genetic effect (σ²a = 20.0, σ²e = 40.0)
a = RandomEffect(Z, sparse(Ai), 20.0; name = "animal")

# Assemble and solve MME
m = mme(wwg, X, [a]; σ²e = 40.0)
sol = solve_mme(m)
b, (u,) = solutions(m, sol)

# Accuracies
pev_animal = pev(m, 1)
r2_animal  = reliability(m, 1)
```

---

## Documentation

To build the static documentation locally:

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs --startup-file=no docs/make.jl
```

The output site will be in `docs/build/`.
