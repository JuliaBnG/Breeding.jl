# API reference

## Package

```@docs
Breeding.Breeding
```

## Mixed model equations (MME)

```@docs
RandomEffect
MME
mme
solve_mme
solutions
lhs_inverse
pev
reliability
```

## Solvers and iteration-on-data (IOD)

```@docs
pcg
gauss_seidel
jacobi
AnimalModelIOD
rhs
```

## Design matrices and multi-trait helpers

```@docs
incidence
legendre
linear_spline
mt_blocks
mt_stack
mt_rinv
associate_incidence
mtpblup
```

## Variance component estimation

```@docs
reml
REMLResult
gibbs
gibbs_mt
gibbs_sweep!
GibbsResult
henderson3
sire_contrasts
fit_mean_squares
covariance_function
```

## Genomic prediction and Bayesian alphabet

```@docs
center_genotypes
snp_from_gblup
gblup_index
base_allele_frequencies
pseudo_snps
bayes_snp
BayesResult
```

## Threshold and categorical models

```@docs
threshold_model
ThresholdResult
category_probabilities
joint_binary_quantitative
JointBinaryResult
```

## Pedigree partitions and reduced animal model (RAM)

```@docs
mendelian_precision
ebv_partition
ebv_partition_mt
ram_design
ram_backsolve
```

## Simulation and reproduction

```@docs
gene_drop
sum_map
sumMap
sampleID
sampleLoci
eqtl
tbv!
phenotype!
```
