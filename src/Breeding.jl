"""
    Breeding

Quantitative-genetics tools for breeding-programme simulation, genetic
evaluation and genomic prediction from pedigree and marker data.

# Features
- **Mixed-model equations (MME)**: Henderson's MME with any number of
  random terms, each with correlated components (``G_0 \\otimes K``):
  multiple traits with missing records, maternal and
  social/associative effects, random regressions on Legendre
  polynomials or linear splines; PEV and reliabilities.
- **Solvers and iteration on data**: sparse Cholesky, preconditioned
  conjugate gradients (`pcg`), Gauss–Seidel and SOR (`gauss_seidel`
  with `ω`), Jacobi (`jacobi`), and matrix-free iteration on data
  (`AnimalModelIOD`).
- **Variance components**: REML (`reml`, average-information or EM
  updates), Gibbs samplers (`gibbs`, `gibbs_mt`), Henderson's method 3
  (`henderson3`), canonical sire contrasts (`sire_contrasts`,
  `fit_mean_squares`) and covariance functions (`covariance_function`).
- **Genomic prediction**: genotype centring (`center_genotypes`), GBLUP
  in selection-index form (`gblup_index`), SNP effects from GBLUP
  (`snp_from_gblup`), base allele frequencies (`base_allele_frequencies`),
  pseudo-SNPs from phased haplotypes (`pseudo_snps`), and Bayes A, B, C
  and Cπ (`bayes_snp`).
- **Threshold models**: probit threshold models solved by Fisher
  scoring (`threshold_model`), category probabilities
  (`category_probabilities`), and joint analysis of a quantitative and
  a binary trait (`joint_binary_quantitative`).
- **Partitions and reduced animal model**: Mendelian-sampling
  precisions (`mendelian_precision`), PA/YD/PC and DYD partitions of
  EBV (`ebv_partition`, `ebv_partition_mt`), and the reduced animal
  model (`ram_design`, `ram_backsolve`).
- **Simulation**: gene dropping with Poisson-distributed crossovers
  (`gene_drop`), sampling of individuals and loci (`sampleID`,
  `sampleLoci`), QTL effects (`eqtl`), true breeding values (`tbv!`) and
  phenotypes (`phenotype!`).
"""
module Breeding

using BnGStructs
using FisherWright
using RelationshipMatrices
using DataFrames
using Distributions
using LinearAlgebra
using Statistics
using Random
using SparseArrays
using StatsBase

include("util.jl")
include("reproduction/reproduce.jl")
include("prediction/predict.jl")
include("vc/vc.jl")
include("trait/trait.jl")
include("ocs/ocs.jl")

export sampleID, sampleLoci, eqtl, tbv!, phenotype!, mtpblup
export gene_drop, sum_map, sumMap
export incidence, legendre, linear_spline, mt_blocks, mt_stack, mt_rinv, associate_incidence
export RandomEffect, MME, mme, solve_mme, solutions, lhs_inverse, pev, reliability
export jacobi, gauss_seidel, pcg, AnimalModelIOD, rhs
export covariance_function, henderson3, sire_contrasts, fit_mean_squares, reml, REMLResult
export gibbs_sweep!, gibbs, gibbs_mt, GibbsResult
export center_genotypes, snp_from_gblup, gblup_index, base_allele_frequencies, pseudo_snps
export bayes_snp, BayesResult
export threshold_model, ThresholdResult, category_probabilities
export joint_binary_quantitative, JointBinaryResult
export mendelian_precision, ebv_partition, ebv_partition_mt, ram_design, ram_backsolve

end # module Breeding
