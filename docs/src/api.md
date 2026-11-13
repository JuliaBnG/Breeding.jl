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
pev_candidates
selected_inverse
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
nested
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
ResidualStructure
Homogeneous
HetClasses
MultiTrait
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
zpz!
zpz_hap!
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

## Optimum contribution selection (OCS & ROCS)

```@docs
TM1997
OC_fixsex2
rocs
ROCSResult
```
