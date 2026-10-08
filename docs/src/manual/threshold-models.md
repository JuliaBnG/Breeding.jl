# Threshold & Categorical Models

Many economically critical agricultural traits—such as calving ease, disease incidence, survival,
and fertility—are observed as binary outcomes or ordered categorical scores rather than continuous
measurements.

`Breeding.jl` implements the liability threshold model (Wright, 1934;
Gianola and Foulley, 1983), solved for the posterior mode by Fisher
scoring, and the joint analysis of a continuous and a binary trait.

---

## Ordered Categorical Threshold Model

Following Gianola and Foulley (1983) and Harville and Mee (1984), an observed categorical response
$Y \in \{1, 2, \dots, m\}$ is generated when an unobserved continuous liability $\ell$ falls between
fixed thresholds:

```math
Y = k \iff t_{k-1} < \ell \le t_k, \quad \text{with } t_0 = -\infty, \ t_m = +\infty
```

Under the probit link function, the conditional probability of category $k$ given liability mean
$\eta = \mathbf{x}'\mathbf{b} + \mathbf{z}'\mathbf{u}$ is:

```math
P(Y = k \mid \eta) = \Phi(t_k - \eta) - \Phi(t_{k-1} - \eta)
```

where $\Phi(\cdot)$ is the standard normal cumulative distribution function. The residual variance
of the liability is set to $1$ for identification.

### Newton–Raphson Fisher Scoring

[`threshold_model`](@ref) iteratively solves the augmented scoring equations:

```math
\begin{bmatrix}
\mathbf{Q} & \mathbf{L}'\mathbf{X} & \mathbf{L}'\mathbf{Z} \\
\mathbf{X}'\mathbf{L} & \mathbf{X}'\mathbf{W}\mathbf{X} & \mathbf{X}'\mathbf{W}\mathbf{Z} \\
\mathbf{Z}'\mathbf{L} & \mathbf{Z}'\mathbf{W}\mathbf{X} & \mathbf{Z}'\mathbf{W}\mathbf{Z} + \lambda \mathbf{K}^{-1}
\end{bmatrix}
\begin{bmatrix}
\Delta \mathbf{t} \\ \Delta \mathbf{b} \\ \Delta \mathbf{u}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{p} \\ \mathbf{X}'\mathbf{v} \\ \mathbf{Z}'\mathbf{v} - \lambda \mathbf{K}^{-1}\mathbf{u}
\end{bmatrix}
```

```julia
# N: contingency table of subclass counts across m categories
res = threshold_model(N, X, Z, ainv(ped), λ; t0 = [0.5, 1.2], zero = [1])

# Extract threshold and breeding value solutions:
thresholds = res.t
ebv_liability = res.u
```

### Computing Predicted Probabilities

The probabilities of falling into each category for individuals with liability means $\boldsymbol{\eta}$
are computed via [`category_probabilities`](@ref):

```julia
P = category_probabilities(res.t, X * res.b + Z * res.u)
```

---

## Joint Continuous and Binary Analysis

In breeding programs, continuous performance traits (e.g., birth weight) often exhibit genetic and
environmental correlations with binary threshold traits (e.g., calving difficulty).

[`joint_binary_quantitative`](@ref) implements Foulley et al.'s (1983) joint mixed model (Mrode Section 15.3).
The liability of the binary trait is conditioned on the observed continuous record via the residual
regression:

```math
\beta = \frac{r_{12}}{\sigma_{e1}\sqrt{1 - r_{12}^2}}
```

This decouples the residual errors, enabling stable iterative estimation of fixed and genetic effects
across both scales.

```julia
# d_binary = 1 for the category with probability Φ(μ), e.g. difficult calving
res_joint = joint_binary_quantitative(y_quant, d_binary, X1, Z1, X2, Z2, ainv(ped), G_cov, σ²e1, r12)
```
