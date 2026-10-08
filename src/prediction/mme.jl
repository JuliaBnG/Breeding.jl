"""
    RandomEffect(Z, Kinv, σ²; name = "")
    RandomEffect(Zs::Vector, Kinv, G0; name = "")

A random term of a linear mixed model with ``m`` correlated components,
each with ``q`` levels:

```math
\\mathbf{u} = [\\mathbf{u}_1' \\cdots \\mathbf{u}_m']', \\quad
\\operatorname{var}(\\mathbf{u}) = \\mathbf{G}_0 \\otimes \\mathbf{K},
```

contributing ``\\sum_c \\mathbf{Z}_c \\mathbf{u}_c`` to the records.
`Zs[c]` (`nobs × q`) relates records to component `c`, `Kinv` is
``\\mathbf{K}^{-1}`` (`q × q`; e.g. ``A^{-1}``, ``G^{-1}``, ``H^{-1}``,
``D^{-1}`` or `I`) and `G0` is the `m × m` covariance among components.

This covers, among others:

- one component: animal, sire, permanent environment, litter, SNP, …;
- direct + maternal (or direct + social) genetic effects, where `Zs`
  are the calf and dam (or group-mate) incidence matrices;
- multiple traits, where `Zs[t]` maps trait `t` records (see
  [`mt_blocks`](@ref));
- random regressions, where `Zs[c] = Diagonal(Φ[:, c]) * Z`.

Solutions are ordered component-major (all levels of component 1,
then component 2, …), matching ``\\mathbf{G}_0^{-1} \\otimes
\\mathbf{K}^{-1}``.
"""
struct RandomEffect{TK<:AbstractMatrix{Float64}}
    Z::SparseMatrixCSC{Float64,Int}  # nobs × (m⋅q), component-major
    Kinv::TK                         # q × q
    G0::Matrix{Float64}              # m × m
    name::String
end

_kinv(Kinv::UniformScaling, q) = sparse(Float64(Kinv.λ) * I, q, q)
_kinv(Kinv::SparseMatrixCSC, q) = SparseMatrixCSC{Float64,Int}(Kinv)
_kinv(Kinv::AbstractMatrix, q) = Matrix{Float64}(Kinv)

function RandomEffect(Zs::AbstractVector{<:AbstractMatrix}, Kinv, G0::AbstractMatrix;
                      name = "")
    m = length(Zs)
    size(G0) == (m, m) || throw(DimensionMismatch("G0 must be $m × $m"))
    q = size(Zs[1], 2)
    all(size(Z, 2) == q for Z in Zs) ||
        throw(DimensionMismatch("all components need the same number of levels"))
    K = _kinv(Kinv, q)
    size(K) == (q, q) || throw(DimensionMismatch("Kinv must be $q × $q"))
    Z = SparseMatrixCSC{Float64,Int}(reduce(hcat, sparse.(Zs)))
    RandomEffect(Z, K, Matrix{Float64}(G0), string(name))
end

RandomEffect(Z::AbstractMatrix, Kinv, σ²::Real; name = "") =
    RandomEffect([Z], Kinv, fill(Float64(σ²), 1, 1); name = name)

nlevels(r::RandomEffect) = size(r.Kinv, 1)
ncomponents(r::RandomEffect) = size(r.G0, 1)

"""
    MME

Mixed-model equations ``\\mathbf{C}\\hat{\\theta} = \\mathbf{r}`` built
by [`mme`](@ref). Fields: `lhs`, `rhs`, the equation ranges `blocks`
(fixed effects first, then one range per random effect), `effects`, and
`scale` – the residual variance the equations were multiplied by (`1`
when a full ``R^{-1}`` was supplied), so that
``\\operatorname{PEV} = C^{-1} \\cdot`` `scale`.
"""
struct MME{TE}
    lhs::SparseMatrixCSC{Float64,Int}
    rhs::Vector{Float64}
    blocks::Vector{UnitRange{Int}}
    effects::TE
    scale::Float64
end

"""
    mme(y, X, effects; σ²e = nothing, Rinv = nothing) -> MME

Set up Henderson's mixed-model equations for

```math
\\mathbf{y} = \\mathbf{Xb} + \\sum_k \\mathbf{Z}_k \\mathbf{u}_k + \\mathbf{e},
\\quad \\operatorname{var}(\\mathbf{u}_k) = \\mathbf{G}_{0k} \\otimes \\mathbf{K}_k,
\\quad \\operatorname{var}(\\mathbf{e}) = \\mathbf{R}.
```

- With `σ²e` (``\\mathbf{R} = \\mathbf{I}\\sigma^2_e``) the equations are
  scaled by ``\\sigma^2_e``, as printed in Mrode:
  ``\\mathbf{W'W} + \\sigma^2_e \\oplus_k (\\mathbf{G}_{0k}^{-1} \\otimes
  \\mathbf{K}_k^{-1})``, giving the familiar variance ratios ``\\alpha``.
- With `Rinv` (a matrix, or a vector of diagonal weights) the
  unscaled ``\\mathbf{W'R^{-1}W} + \\oplus_k (\\mathbf{G}_{0k}^{-1}
  \\otimes \\mathbf{K}_k^{-1})`` is built. Use it for heterogeneous
  residuals, reduced animal models, or multiple traits
  ([`mt_rinv`](@ref)).

`X` may have zero columns; `effects` is a vector of
[`RandomEffect`](@ref)s (possibly empty). All products are sparse.
"""
function mme(y::AbstractVector, X::AbstractMatrix, effects::AbstractVector{<:RandomEffect} = RandomEffect[];
             σ²e = nothing, Rinv = nothing)
    (σ²e === nothing) == (Rinv === nothing) &&
        throw(ArgumentError("give exactly one of σ²e or Rinv"))
    n = length(y)
    size(X, 1) == n || throw(DimensionMismatch("X has $(size(X, 1)) rows, y has $n"))
    for r in effects
        size(r.Z, 1) == n || throw(DimensionMismatch("Z of '$(r.name)' has wrong rows"))
    end
    W = SparseMatrixCSC{Float64,Int}(hcat(sparse(X), (r.Z for r in effects)...))
    yv = Vector{Float64}(y)
    if Rinv === nothing
        scale = Float64(σ²e)
        lhs = SparseMatrixCSC{Float64,Int}(W'W)
        rhs = W'yv
    else
        scale = 1.0
        Ri = Rinv isa AbstractVector ? Diagonal(Vector{Float64}(Rinv)) : Rinv
        WR = SparseMatrixCSC{Float64,Int}(W' * Ri)
        lhs = WR * W
        rhs = WR * yv
    end
    blocks = UnitRange{Int}[1:size(X, 2)]
    pen = SparseMatrixCSC{Float64,Int}[spzeros(size(X, 2), size(X, 2))]
    for r in effects
        start = last(blocks[end]) + 1
        push!(blocks, start:start+size(r.Z, 2)-1)
        Gi = inv(Symmetric(r.G0)) .* scale
        push!(pen, kron(sparse(Gi), r.Kinv isa SparseMatrixCSC ? r.Kinv : sparse(r.Kinv)))
    end
    lhs += blockdiag(pen...)
    MME(dropzeros!(lhs), rhs, blocks, effects, scale)
end

"""
    _keep(n, zero) -> Vector{Int}

Indices of the equations that are not constrained to zero.
"""
_keep(n, zero) = isempty(zero) ? collect(1:n) : setdiff(1:n, zero)

"""
    solve_mme(m::MME; zero = Int[], method = :direct, kw...) -> Vector{Float64}
    solve_mme(lhs, rhs; zero = Int[], method = :direct, kw...)

Solve the mixed-model equations. Equations listed in `zero` are
removed and their solutions set to zero; this imposes the usual
``\\hat{b}_k = 0`` constraints when ``\\mathbf{X}`` is not of full
column rank (e.g. one level of each extra fixed factor, or one unknown
parent group).

`method`:
- `:direct` – sparse Cholesky (CHOLMOD) of the positive definite
  system; falls back to sparse LU if Cholesky fails;
- `:pinv` – dense Moore–Penrose solution (small, rank-deficient
  systems only);
- `:pcg`, `:gauss_seidel`, `:sor`, `:jacobi` – iterative solvers of
  Chapter 19 ([`pcg`](@ref), [`gauss_seidel`](@ref),
  [`jacobi`](@ref)); keywords are passed on.
"""
solve_mme(m::MME; kw...) = solve_mme(m.lhs, m.rhs; kw...)

function solve_mme(lhs::AbstractMatrix, rhs::AbstractVector; zero = Int[],
                   method::Symbol = :direct, kw...)
    n = length(rhs)
    keep = _keep(n, zero)
    A = lhs[keep, keep]
    b = rhs[keep]
    x = zeros(n)
    x[keep] = if method == :direct
        _direct(A, b)
    elseif method == :pinv
        pinv(Matrix(A)) * b
    elseif method == :pcg
        first(pcg(A, b; kw...))
    elseif method == :gauss_seidel
        first(gauss_seidel(A, b; kw...))
    elseif method == :sor
        first(gauss_seidel(A, b; merge((; ω = 1.2), kw)...))
    elseif method == :jacobi
        first(jacobi(A, b; kw...))
    else
        throw(ArgumentError("unknown method $method"))
    end
    x
end

function _direct(A::SparseMatrixCSC, b)
    F = cholesky(Symmetric(A); check = false)
    issuccess(F) && return F \ b
    lu(A) \ b
end
_direct(A::AbstractMatrix, b) = Symmetric(Matrix(A)) \ b

"""
    lhs_inverse(m::MME; zero = Int[]) -> Matrix{Float64}
    lhs_inverse(lhs; zero = Int[])

Dense (generalized) inverse of the MME coefficient matrix, with zero
rows and columns for the constrained equations. Intended for the small
examples where the full ``\\mathbf{C}^{-1}`` is inspected; use
[`pev`](@ref) for prediction-error variances of large models.
"""
lhs_inverse(m::MME; zero = Int[]) = lhs_inverse(m.lhs; zero = zero)

function lhs_inverse(lhs::AbstractMatrix; zero = Int[])
    n = size(lhs, 1)
    keep = _keep(n, zero)
    C = zeros(n, n)
    C[keep, keep] = inv(Symmetric(Matrix(lhs[keep, keep])))
    C
end

"""
    pev(m::MME, k::Integer; zero = Int[]) -> Matrix{Float64}

Prediction-error variances ``\\operatorname{var}(u - \\hat{u})`` of the
`k`th random effect, as a `q × m` matrix (levels × components). They
are the diagonal of the corresponding block of ``\\mathbf{C}^{-1}``
times `m.scale`. The coefficient matrix is factorized once and only the
required columns of the inverse are solved for.
"""
function pev(m::MME, k::Integer; zero = Int[])
    r = m.effects[k]
    rng = m.blocks[k+1]
    n = length(m.rhs)
    keep = _keep(n, zero)
    pos = zeros(Int, n)
    pos[keep] = 1:length(keep)
    A = m.lhs[keep, keep]
    F = cholesky(Symmetric(A); check = false)
    fac = issuccess(F) ? F : lu(A)
    d = zeros(length(rng))
    e = zeros(length(keep))
    for (i, j) in enumerate(rng)
        pos[j] == 0 && continue
        fill!(e, 0)
        e[pos[j]] = 1
        d[i] = (fac\e)[pos[j]]
    end
    reshape(d .* m.scale, nlevels(r), ncomponents(r))
end

"""
    reliability(m::MME, k::Integer; zero = Int[], Kdiag = nothing) -> Matrix{Float64}

Reliability ``r^2 = 1 - \\operatorname{PEV}/\\sigma^2_u`` of the
levels of random effect `k`, per component (``\\sigma^2_u`` is the
diagonal of `G0` times the diagonal of ``\\mathbf{K}``). Without
`Kdiag`, ``\\operatorname{diag}(\\mathbf{K})`` is obtained by densely
inverting ``\\mathbf{K}^{-1}``, which is only practical for small
models; for ``\\mathbf{K} = A`` pass `Kdiag = nrm_diag(ped)` (i.e.
``1 + F_i``). For non-inbred animals this is Henderson's
``1 - d_i\\alpha``.
"""
function reliability(m::MME, k::Integer; zero = Int[], Kdiag = nothing)
    r = m.effects[k]
    P = pev(m, k; zero = zero)
    kd = Kdiag === nothing ? diag(inv(Matrix(r.Kinv))) : Kdiag
    1 .- P ./ (kd * diag(r.G0)')
end

"""
    solutions(m::MME, x) -> (fixed = b, random = [u₁, u₂, …])

Split a solution vector into fixed effects and, per random effect, a
`q × m` matrix (levels in rows, components in columns).
"""
function solutions(m::MME, x::AbstractVector)
    fixed = x[m.blocks[1]]
    random = [reshape(x[m.blocks[k+1]], nlevels(r), ncomponents(r))
              for (k, r) in enumerate(m.effects)]
    (fixed = fixed, random = random)
end
