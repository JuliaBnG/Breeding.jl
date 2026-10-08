# Linear Mixed Models & BLUP

Linear mixed models form the bedrock of quantitative genetic evaluation, allowing simultaneous
estimation of fixed environmental effects (Best Linear Unbiased Estimation, BLUE) and prediction
of random additive genetic effects (Best Linear Unbiased Prediction, BLUP).

`Breeding.jl` provides an extensible, sparse-matrix formulation of Henderson's Mixed Model Equations
(MME) capable of accommodating single-trait animal models, maternal effects, random regression,
multi-trait models with arbitrary patterns of missing phenotypes, and reduced animal models.

---

## Henderson's Mixed Model Equations

Consider the general linear mixed model:

```math
\mathbf{y} = \mathbf{Xb} + \sum_{k} \mathbf{Z}_k \mathbf{u}_k + \mathbf{e}
```

where:
- $\mathbf{y}$ is the $n \times 1$ vector of observed phenotypic records;
- $\mathbf{b}$ is the $p \times 1$ vector of fixed effects with design matrix $\mathbf{X}$;
- $\mathbf{u}_k$ is the random vector for effect $k$ with $m_k$ components and $q_k$ levels per component;
- $\operatorname{var}(\mathbf{u}_k) = \mathbf{G}_{0k} \otimes \mathbf{K}_k$, where $\mathbf{G}_{0k}$ is the $m_k \times m_k$ covariance among components and $\mathbf{K}_k$ is the $q_k \times q_k$ relationship matrix (e.g., pedigree numerator relationship $\mathbf{A}$, genomic relationship $\mathbf{G}$, or single-step $\mathbf{H}$);
- $\mathbf{e}$ is the residual vector with $\operatorname{var}(\mathbf{e}) = \mathbf{R}$.

Henderson's Mixed Model Equations take the form:

```math
\begin{bmatrix}
\mathbf{X}'\mathbf{R}^{-1}\mathbf{X} & \mathbf{X}'\mathbf{R}^{-1}\mathbf{Z} \\
\mathbf{Z}'\mathbf{R}^{-1}\mathbf{X} & \mathbf{Z}'\mathbf{R}^{-1}\mathbf{Z} + \bigoplus_k (\mathbf{G}_{0k}^{-1} \otimes \mathbf{K}_k^{-1})
\end{bmatrix}
\begin{bmatrix}
\hat{\mathbf{b}} \\ \hat{\mathbf{u}}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{X}'\mathbf{R}^{-1}\mathbf{y} \\
\mathbf{Z}'\mathbf{R}^{-1}\mathbf{y}
\end{bmatrix}
```

When residuals are homogeneous and independent ($\mathbf{R} = \mathbf{I}\sigma^2_e$), the equations can be scaled by $\sigma^2_e$, yielding the familiar variance ratio $\alpha = \sigma^2_e / \sigma^2_a$:

```math
\begin{bmatrix}
\mathbf{X}'\mathbf{X} & \mathbf{X}'\mathbf{Z} \\
\mathbf{Z}'\mathbf{X} & \mathbf{Z}'\mathbf{Z} + \alpha \mathbf{K}^{-1}
\end{bmatrix}
\begin{bmatrix}
\hat{\mathbf{b}} \\ \hat{\mathbf{u}}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{X}'\mathbf{y} \\
\mathbf{Z}'\mathbf{y}
\end{bmatrix}
```

---

## Defining Random Effects: `RandomEffect`

A random term is represented by `RandomEffect`. It handles:
1. **Single-component terms**: Animal additive genetic effect, sire effect, permanent environment, or litter effect.
2. **Multi-component terms**:
   - Direct and maternal genetic effects: $\mathbf{Z} = [\mathbf{Z}_{\text{direct}} \quad \mathbf{Z}_{\text{maternal}}]$ with covariance $\mathbf{G}_0 = \begin{bmatrix} \sigma^2_a & \sigma_{am} \\ \sigma_{am} & \sigma^2_m \end{bmatrix}$.
   - Social / indirect genetic effects (Bijma et al., 2007; Muir, 2005).
   - Multiple traits stacked trait-major.
   - Random regression models with orthogonal Legendre polynomials or linear splines.

```julia
# Single-component random effect
a = RandomEffect(Z, ainv(ped), σ²a; name = "animal")

# Direct-maternal model with covariance G0
G0 = [20.0 -4.0; -4.0 10.0]
r = RandomEffect([Z_direct, Z_maternal], ainv(ped), G0; name = "direct_maternal")
```

---

## Design Matrix Construction

`Breeding.jl` provides specialized functions for setting up incidence matrices:

- [`incidence`](@ref): Build sparse design matrices from factor vectors or index arrays.
- [`legendre`](@ref): Evaluate normalized Legendre polynomials up to order $k-1$ standardized to $[-1, 1]$ (Kirkpatrick et al., 1990; Mrode Appendix G).
- [`linear_spline`](@ref): Construct linear spline basis functions across knot intervals (Misztal, 2006).
- [`associate_incidence`](@ref): Set up associative/social genetic effect design matrices relating each animal to its pen- or cage-mates with optional dilution weighting.
- [`mt_stack`](@ref), [`mt_blocks`](@ref), [`mt_rinv`](@ref): Handle multivariate phenotypes with missing records, computing exact block-inverse residual matrices $\mathbf{R}^{-1}$.

---

## Setting Up and Solving MME

Assemble equations using [`mme`](@ref) and solve with [`solve_mme`](@ref):

```julia
# Homogeneous residual variance:
m = mme(y, X, [a]; σ²e = 40.0)

# Heterogeneous or multivariate residual matrix:
m = mme(y, X, [a]; Rinv = R_inv)

# Solve using sparse Cholesky factorization:
sol = solve_mme(m)

# Separate fixed and random solutions:
b, (u,) = solutions(m, sol)
```

Fixed-effect rank deficiencies (e.g. overparameterized contemporary groups) can be handled by constraining chosen equations to zero via `zero = [...]`:

```julia
sol = solve_mme(m; zero = [1])
```

---

## Prediction Error Variance (PEV) and Reliability

The prediction error variance of animal $i$ is:

```math
\operatorname{PEV}_i = \operatorname{var}(u_i - \hat{u}_i) = C^{ii}
```

where $C^{ii}$ is the diagonal element of the inverse of the unscaled
coefficient matrix ($\mathbf{R}^{-1}$ form) corresponding to animal
$i$; for equations scaled by $\sigma^2_e$ it is $C^{ii}\sigma^2_e$,
which `pev` applies automatically. The theoretical reliability is:

```math
r^2_i = 1 - \frac{\operatorname{PEV}_i}{(1 + F_i)\sigma^2_a}
```

[`pev`](@ref) factorizes the coefficient matrix once and solves only
for the required columns of its inverse. [`reliability`](@ref) also
needs $\operatorname{diag}(\mathbf{K})$; pass it as `Kdiag` (e.g.
`nrm_diag(ped)` for $1 + F_i$), otherwise it is obtained by densely
inverting $\mathbf{K}^{-1}$, which suits small models only:

```julia
# PEV and reliability for random effect 1:
P = pev(m, 1)
r2 = reliability(m, 1)
```

---

## Reduced Animal Model (RAM)

In large populations with many unrecorded or non-parent animals, Quaas and Pollak's (1980)
Reduced Animal Model sets up equations solely for parents. Non-parent records are mapped
to their parents with variance augmented by Mendelian sampling factors:

```julia
# Generate RAM design and inverse residual weights
W, rinv, parents, Api = ram_design(ped, records, σ²a, σ²e)

# Solve parent equations
m_ram = mme(y, X, [RandomEffect(W, Api, σ²a)]; Rinv = rinv)
sol_ram = solve_mme(m_ram)
b_ram, (ap,) = solutions(m_ram, sol_ram)

# Back-solve non-parent EBVs from records and parent averages
a_all = ram_backsolve(ped, records, y - X * b_ram, vec(ap), parents, σ²a, σ²e)
```

---

## EBV Partitioning

VanRaden and Wiggans' (1991) breeding value decomposition partitions an animal's EBV
into parent average ($\mathrm{PA}$), yield deviation ($\mathrm{YD}$), and progeny contribution ($\mathrm{PC}$):

```math
\hat{a}_i = n_1 \mathrm{PA}_i + n_2 \mathrm{YD}_i + n_3 \mathrm{PC}_i = w_1 \mathrm{PA}_i + w_2 \mathrm{YD}_i + w_3 \mathrm{DYD}_i
```

Both univariate ([`ebv_partition`](@ref)) and multi-trait ([`ebv_partition_mt`](@ref)) formulations
are implemented.
