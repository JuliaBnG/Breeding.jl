# Iterative Solvers & Iteration-on-Data

When genetic evaluations scale to hundreds of thousands or millions of animals, forming and factoring
the mixed model coefficient matrix $\mathbf{C}$ becomes computationally infeasible due to fill-in
and memory constraints.

`Breeding.jl` provides standard stationary iterative solvers, preconditioned Krylov subspace solvers,
and matrix-free **Iteration-on-Data** (IOD) operators that read pedigree and record streams directly.

---

## Stationary Iterative Solvers

### Jacobi Iteration

The Jacobi method updates each equation using the solutions from the previous iteration:

```math
\mathbf{x}^{(r+1)} = \omega \mathbf{M}^{-1}(\mathbf{b} - \mathbf{A}\mathbf{x}^{(r)}) + \mathbf{x}^{(r)}
```

where $\mathbf{M} = \operatorname{diag}(\mathbf{A})$. For a symmetric
positive definite $\mathbf{A}$, damped Jacobi converges when
$0 < \omega < 2/\rho(\mathbf{M}^{-1}\mathbf{A})$; plain Jacobi
($\omega = 1$) can diverge on animal-model equations, which is why
Mrode uses $\omega = 0.8$ for the animal equations (`ω` may be a vector
of per-equation factors). Second-order Jacobi adds a momentum term
$\omega(\mathbf{x}^{(r)} - \mathbf{x}^{(r-1)})$.

```julia
x, iters, conv = jacobi(A, b; ω = 0.8, tol = 1e-8)
```

### Gauss–Seidel and SOR

Gauss–Seidel immediately incorporates updated solutions into subsequent equations within the same round. For symmetric sparse systems stored in compressed sparse column (CSC) format, row $i$ is accessed directly from column $i$:

```julia
x, iters, conv = gauss_seidel(A, b; tol = 1e-9)

# Successive Over-Relaxation (SOR with ω = 1.2):
x, iters, conv = gauss_seidel(A, b; ω = 1.2, tol = 1e-9)
```

---

## Preconditioned Conjugate Gradient (PCG)

For symmetric positive definite mixed-model systems, the
preconditioned conjugate gradient (PCG; Lidauer et al., 1999; Mrode
Section 19.6) usually converges in fewer rounds than stationary methods:

```julia
# PCG on sparse matrix with diagonal preconditioner:
x, iters, conv = pcg(A, b; M = diag(A), tol = 1e-12)
```

The operator `A` in [`pcg`](@ref) can be an explicit matrix or an in-place linear mapping function `A(v, d)` that evaluates $\mathbf{v} \leftarrow \mathbf{A}\mathbf{d}$.

---

## Iteration-on-Data (IOD)

The Iteration-on-Data algorithm (Schaeffer and Kennedy, 1986) evaluates mixed model equations
without ever forming or storing non-zero elements of $\mathbf{C}$.

### `AnimalModelIOD`

[`AnimalModelIOD`](@ref) stores only:
1. Phenotypic records and fixed factor levels;
2. Animal indices, sire and dam pointers;
3. Mendelian sampling precision factors $\delta_i = \alpha / (1 - \tfrac14 \sum_{\text{known } p} (1 + F_p))$;
4. Progeny and mate lists (CSR format).

It supports:
- any number of fixed factors;
- inbreeding coefficients $F_i$ (computed from the pedigree by the
  `DataFrame` constructor);
- unknown parent groups (phantom parents; Westell et al., 1988),
  appended after the animal equations. Groups are coded `N + g` in the
  `sire`/`dam` vectors of the vector constructor; the `DataFrame`
  constructor does not handle groups.

```julia
# Construct IOD model from DataFrame pedigree
iod = AnimalModelIOD(y, [fixed_factor_1, fixed_factor_2], animal_idx, ped, α)

# Compute right-hand side directly:
r = rhs(iod)

# Extract diagonal preconditioner:
M = diag(iod)

# Solve via matrix-free PCG:
x, iters, conv = pcg(iod, r; M = M, tol = 1e-10)

# Or solve via Gauss-Seidel iteration on data:
x_gs, iters, conv = gauss_seidel(iod; tol = 1e-10)
```
