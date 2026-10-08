# Variance Component Estimation

Reliable estimation of genetic and environmental variance components is essential for accurate
genetic evaluation, breeding program design, and predicting selection response.

`Breeding.jl` provides modern likelihood-based methods (REML), fully Bayesian Markov chain Monte
Carlo samplers (Gibbs), classical ANOVA-based estimation (Henderson's Method 3), and longitudinal
covariance functions.

---

## Restricted Maximum Likelihood (REML)

[`reml`](@ref) implements Restricted Maximum Likelihood (Patterson and Thompson, 1971) for general
mixed models:

```math
\mathbf{y} = \mathbf{Xb} + \sum_{k} \mathbf{Z}_k \mathbf{u}_k + \mathbf{e}, \quad \operatorname{var}(\mathbf{u}_k) = \mathbf{G}_{0k} \otimes \mathbf{K}_k, \quad \operatorname{var}(\mathbf{e}) = \mathbf{I}\sigma^2_e
```

### Optimization Algorithms

1. **Average Information (AI-REML)**:
   Gilmour et al.'s (1995) algorithm replaces the expected information with the average of observed and expected information, computed efficiently using working variates $\mathbf{f}_i = (\partial \mathbf{V} / \partial \theta_i)\mathbf{P}\mathbf{y}$. Positive-definiteness is ensured via step halving.
2. **EM-type updates (Mrode Eqns 17.5–17.6)**:
   $\sigma^2_e = \mathbf{e}'\mathbf{y} / (n - p)$ and
   $\mathbf{G}_{0k} = (\hat{\mathbf{U}}_k'\mathbf{K}_k^{-1}\hat{\mathbf{U}}_k + \mathbf{T}_k)/q_k$.
   Estimates stay positive definite, but convergence is slow: in Mrode's
   Section 17.7 AI-REML converges in 5 rounds, while these updates are
   still changing in the fourth decimal after 1000.

```julia
# Fit AI-REML:
res = reml(y, X, [RandomEffect(Z, ainv(ped), 0.2)]; σ²e = 0.4, method = :ai)

# Estimates and standard errors:
σ²e_est = res.σ²e
G0_est  = res.G0[1]
vc_cov  = res.AIinv  # approximate sampling covariance (inverse AI matrix)
```

---

## Bayesian Gibbs Sampling

### Univariate Mixed Model (`gibbs`)

[`gibbs`](@ref) performs Gibbs sampling for location parameters and variance components (Wang et al., 1993; Sorensen and Gianola, 2002):

```math
\sigma^2_e \mid \cdot \sim \frac{\mathbf{e}'\mathbf{e} + \nu_e S_e}{\chi^2_{n + \nu_e}}, \quad \mathbf{G}_{0k} \mid \cdot \sim \operatorname{IW}(\nu_{uk} + q_k, \mathbf{U}_k'\mathbf{K}_k^{-1}\mathbf{U}_k + \mathbf{V}_{uk})
```

Location sweeps ([`gibbs_sweep!`](@ref)) are performed directly on the MME system. Fixing variance components (`estimate_variances = false`) allows sampling posterior EBVs whose means converge to BLUP.

```julia
res = gibbs(y, X, [a]; σ²e = 40.0, niter = 10_000, burnin = 1_000, thin = 5)
```

### Multivariate Block Gibbs (`gibbs_mt`)

[`gibbs_mt`](@ref) implements Jensen et al.'s (1994) multi-trait block Gibbs sampler. Incomplete multivariate records are imputed in each round via **data augmentation** conditioned on observed traits and current covariance estimates.

```julia
res_mt = gibbs_mt(Y_pheno, X, Z, ainv(ped); R0 = R_init, G0 = G_init, niter = 5_000, burnin = 500)
```

---

## Henderson's Method 3 (ANOVA)

For historical benchmarking or balanced sire models, [`henderson3`](@ref) evaluates the classical
sums of squares:

```math
\mathbf{y}'\mathbf{X}(\mathbf{X}'\mathbf{X})^-\mathbf{X}'\mathbf{y}, \quad \mathbf{y}'\mathbf{SZ}(\mathbf{Z}'\mathbf{SZ})^-\mathbf{Z}'\mathbf{Sy}, \quad \mathbf{R} = \mathbf{y}'\mathbf{y} - \mathbf{F} - \mathbf{S}_s
```

Canonical sire contrasts and iterated weighted least squares on mean squares are supported via
[`sire_contrasts`](@ref) and [`fit_mean_squares`](@ref).

---

## Covariance Functions

Longitudinal, test-day, and growth trajectories can be modelled as continuous covariance functions
(Kirkpatrick et al., 1990) of Legendre polynomials:

```math
\mathbf{G} \approx \boldsymbol{\Phi} \mathbf{C} \boldsymbol{\Phi}'
```

[`covariance_function`](@ref) gives the exact coefficient matrix for a
full-order fit and, for a reduced order, fits it by generalized least
squares with Kirkpatrick's $\chi^2$ goodness-of-fit statistic.

```julia
# Fit order-2 covariance function to a 3-age covariance matrix:
cf = covariance_function(G_matrix, [90, 160, 240]; k = 2)
```
