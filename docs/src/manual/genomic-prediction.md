# Genomic Prediction & Bayesian Alphabet

Genomic prediction estimates breeding values using dense genome-wide single nucleotide
polymorphism (SNP) markers or sequence variants.

`Breeding.jl` provides comprehensive utilities for genomic prediction, ranging from marker centring
and selection index formulations to linear GBLUP back-solving and Bayesian alphabet Markov chain
Monte Carlo (MCMC) models.

---

## Genotype Centring & Scaling

Given an $n \times m$ matrix of marker dosages $\mathbf{M}$ ($0, 1, 2$ copies of the counted allele),
VanRaden's (2008) method centres the markers according to allele frequencies:

```math
\mathbf{Z} = \mathbf{M} - 2\mathbf{1}\mathbf{p}'
```

The scaling denominator relates the SNP variance to total additive genetic variance:

```math
k = 2 \sum_{j=1}^{m} p_j(1 - p_j), \quad \sigma^2_g = \frac{\sigma^2_a}{k}
```

```julia
Zc, p, k = center_genotypes(M)
```

---

## Genomic Selection Index (GBLUP)

For reference populations with pre-corrected phenotypes $\mathbf{y}_c$ and variance ratio $\lambda = \sigma^2_e / \sigma^2_a$, [`gblup_index`](@ref) solves the selection index directly for all animals present in $\mathbf{G}$:

```math
\hat{\mathbf{a}} = \mathbf{G}_{\cdot r} (\mathbf{G}_{rr} + \lambda \mathbf{R})^{-1} \mathbf{y}_c
```

Theoretical reliabilities are $b_{ii}/g_{ii}$, with
$\mathbf{B} = \mathbf{G}_{\cdot r}(\mathbf{G}_{rr} + \lambda\mathbf{R})^{-1}\mathbf{G}_{r\cdot}$:

```julia
dgv, rel = gblup_index(G, ref_indices, yc, λ)
```

---

## Back-Solving SNP Effects from GBLUP

Following Strandén and Garrick (2009), SNP effects $\hat{\mathbf{g}}$ can be back-solved from GBLUP breeding values $\hat{\mathbf{a}}$:

```math
\hat{\mathbf{g}} = \frac{1}{k} \mathbf{Z}' \mathbf{G}^{-1} \hat{\mathbf{a}}
```

```julia
g_hat = snp_from_gblup(Zc, inv(G), a_hat; k = k)
```

---

## Base Allele Frequencies (Gengler's Method)

In populations where only a subset of individuals are genotyped, estimating allele frequencies of
the ungenotyped base generation is required to ensure compatibility between pedigree and genomic
relationship matrices.

[`base_allele_frequencies`](@ref) implements Gengler et al.'s (2007) linear mixed model approach,
treating each locus as a continuous trait with pedigree relationship $\mathbf{A}$. The common coefficient
matrix is factorized once, solving all $m$ markers simultaneously:

```julia
p0, μ, d_blup = base_allele_frequencies(M_genotyped, ainv(ped), genotyped_indices; λ = 0.01)
```

---

## Pseudo-SNPs from Phased Haplotypes

When genomic models utilize haplotype alleles across sliding or linkage-disequilibrium blocks
(e.g., Teissier et al., 2020), [`pseudo_snps`](@ref) converts phased haplotypes into pseudo-SNP counts:

```julia
P, alleles = pseudo_snps(haplotypes, [1:10, 11:25, 26:50])
```

---

## Bayesian Alphabet (`bayes_snp`)

`Breeding.jl` implements Gibbs samplers with residual updating (Legarra and Misztal, 2008)
for the regression model $\mathbf{y} = \mathbf{Xb} + \mathbf{Zg} + \mathbf{e}$:

- **BayesA**: Marker-specific variances sampled from inverted chi-squared distributions (Meuwissen et al., 2001).
- **BayesB**: a marker variance of zero with probability $\pi$, else as BayesA; each marker variance is sampled jointly with its effect by Metropolis–Hastings, using the prior as proposal (Meuwissen et al., 2001).
- **BayesC**: Common variance across active markers with indicator variable sampling (Habier et al., 2011).
- **BayesCπ**: BayesC with posterior sampling of the sparsity parameter $\pi \sim \operatorname{Beta}$.

```julia
# Run BayesB with 10,000 iterations and 3,000 burn-in:
res = bayes_snp(y, X, M; method = :B, π = 0.95, niter = 10_000, burnin = 3_000)

# Posterior mean of SNP effects and inclusion probabilities:
g_mean = res.g
pip = res.inclusion
```
