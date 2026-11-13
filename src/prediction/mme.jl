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

function RandomEffect(
    Zs::AbstractVector{<:AbstractMatrix},
    Kinv,
    G0::AbstractMatrix;
    name = "",
)
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
function mme(
    y::AbstractVector,
    X::AbstractMatrix,
    effects::AbstractVector{<:RandomEffect} = RandomEffect[];
    σ²e = nothing,
    Rinv = nothing,
)
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
        push!(blocks, start:(start+size(r.Z, 2)-1))
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
    solve_mme(m::MME; zero = Int[], method = :direct, x0 = nothing, kw...) -> Vector{Float64}
    solve_mme(lhs, rhs; zero = Int[], method = :direct, x0 = nothing, kw...)

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

`x0` (length `length(rhs)`, e.g. the previous round's solutions) is the
starting value of the iterative solvers; entries of `zero` are dropped
from it. The direct methods ignore it.
"""
solve_mme(m::MME; kw...) = solve_mme(m.lhs, m.rhs; kw...)

function solve_mme(
    lhs::AbstractMatrix,
    rhs::AbstractVector;
    zero = Int[],
    method::Symbol = :direct,
    x0 = nothing,
    kw...,
)
    n = length(rhs)
    keep = _keep(n, zero)
    A = lhs[keep, keep]
    b = rhs[keep]
    x = zeros(n)
    if x0 !== nothing
        length(x0) == n || throw(DimensionMismatch("x0 must have length $n"))
        kw = (; kw..., x0 = Vector{Float64}(x0[keep]))
    end
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
    d = _inverse_diagonal(m.lhs; zero = zero)[rng]
    reshape(d .* m.scale, nlevels(r), ncomponents(r))
end

# Diagonal of the (generalized) inverse of a sparse symmetric matrix, with
# zeros for the equations in `zero`. A positive definite matrix goes through
# the selected inverse of its Cholesky factor, whose cost is that of the
# factorisation; otherwise one column is solved per equation with an LU.
function _inverse_diagonal(A::SparseMatrixCSC; zero = Int[])
    n = size(A, 1)
    keep = _keep(n, zero)
    Ak = isempty(zero) ? A : A[keep, keep]
    d = zeros(n)
    F = cholesky(Symmetric(Ak); check = false)
    if issuccess(F)
        Z, _ = selinv(F; depermute = true)
        d[keep] = diag(Z)
    else
        fac = lu(Ak)
        e = zeros(length(keep))
        for i in eachindex(keep)
            e[i] = 1
            d[keep[i]] = (fac\e)[i]
            e[i] = 0
        end
    end
    d
end
_inverse_diagonal(A::AbstractMatrix; zero = Int[]) = diag(lhs_inverse(A; zero = zero))

"""
    selected_inverse(m::MME; zero = Int[]) -> SparseMatrixCSC{Float64,Int}
    selected_inverse(lhs::SparseMatrixCSC; zero = Int[])

Selected elements of the inverse of the MME coefficient matrix: every
element of ``\\mathbf{C}^{-1}`` on the sparsity pattern of the Cholesky
factor of ``\\mathbf{C}``, which contains the pattern of ``\\mathbf{C}``
itself (Takahashi equations, through SelectedInversion.jl). Elements
outside that pattern are not computed and read as zero; they are **not**
zeros of the inverse. Multiply by `m.scale` for prediction error
(co)variances.

This is what large models can afford in place of [`lhs_inverse`](@ref):
memory and time are those of the sparse factorisation. It gives the
diagonal (hence [`pev`](@ref) and [`reliability`](@ref)), the
covariances between an animal and its parents and mates, and the traces
that REML needs. For a dense block among arbitrary candidates use
[`pev_candidates`](@ref).

The coefficient matrix must be positive definite after the equations in
`zero` are removed.
"""
selected_inverse(m::MME; zero = Int[]) = selected_inverse(m.lhs; zero = zero)

function selected_inverse(lhs::SparseMatrixCSC; zero = Int[])
    n = size(lhs, 1)
    keep = _keep(n, zero)
    A = isempty(zero) ? lhs : lhs[keep, keep]
    F = cholesky(Symmetric(A))
    Z, _ = selinv(F; depermute = true)
    Zs = SparseMatrixCSC{Float64,Int}(Z)
    isempty(zero) && return Zs
    I_, J_, V_ = findnz(Zs)
    sparse(keep[I_], keep[J_], V_, n, n)
end

"""
    pev_candidates(m::MME, k, candidates; component = 1, zero = Int[]) -> Matrix{Float64}

Prediction error (co)variance matrix
``\\mathrm{Var}(\\mathbf{u}_c - \\hat{\\mathbf{u}}_c)`` of the levels
`candidates` of random effect `k`, component `component`: the block of
the inverse coefficient matrix for those equations, times the residual
scale. Its diagonal equals the matching elements of [`pev`](@ref). Only
`length(candidates)` right-hand sides are solved, so the full inverse is
never formed.

The result is the matrix `Ω` that [`rocs`](@ref) expects. Equations in
`zero` are constrained to zero as in [`solve_mme`](@ref); a candidate
may not be one of them.
"""
function pev_candidates(
    m::MME,
    k::Integer,
    candidates::AbstractVector{<:Integer};
    component::Integer = 1,
    zero = Int[],
)
    r = m.effects[k]
    q = nlevels(r)
    1 <= component <= ncomponents(r) ||
        throw(ArgumentError("component must be in 1:$(ncomponents(r))"))
    all(c -> 1 <= c <= q, candidates) ||
        throw(ArgumentError("candidates must be levels in 1:$q"))
    offset = first(m.blocks[k+1]) - 1 + (component - 1) * q
    n = length(m.rhs)
    keep = _keep(n, zero)
    pos = zeros(Int, n)
    pos[keep] = 1:length(keep)
    rows = [pos[offset+c] for c in candidates]
    all(>(0), rows) ||
        throw(ArgumentError("a candidate equation is constrained to zero"))
    A = m.lhs[keep, keep]
    F = cholesky(Symmetric(A); check = false)
    fac = issuccess(F) ? F : lu(A)
    E = zeros(length(keep), length(rows))
    for (j, i) in enumerate(rows)
        E[i, j] = 1
    end
    Γ = (fac\E)[rows, :] .* m.scale
    (Γ .+ Γ') ./ 2
end

"""
    pev_candidates(K, ref, candidates, X, λ; σ²a = 1.0, w = nothing) -> Matrix{Float64}

The same matrix from the selection-index (V) form, without inverting
`K`. With records on the animals `ref`, relationship matrix `K` over all
animals, fixed-effect design `X` (one row per record, full column rank)
and ``\\lambda = \\sigma^2_e / \\sigma^2_a``,

```math
\\boldsymbol\\Gamma = \\sigma^2_a \\left(\\mathbf{K}_{cc} -
\\mathbf{K}_{cr}\\mathbf{P}\\mathbf{K}_{rc}\\right), \\qquad
\\mathbf{P} = \\mathbf{V}^{-1} -
\\mathbf{V}^{-1}\\mathbf{X}(\\mathbf{X}'\\mathbf{V}^{-1}\\mathbf{X})^{-1}\\mathbf{X}'\\mathbf{V}^{-1},
```

where ``\\mathbf{V} = \\mathbf{K}_{rr} + \\lambda\\mathbf{I}``, or
``\\mathbf{K}_{rr} + \\lambda\\,\\mathrm{diag}(1/w)`` with record weights
`w` as in [`gblup_index`](@ref). `K` may be `A`, `G`, `H` or an IBD
matrix; it need not be invertible. Candidates may be among `ref`.

This equals the mixed-model form when `K` is invertible and the fixed
effects are the same.
"""
function pev_candidates(
    K::AbstractMatrix,
    ref::AbstractVector{<:Integer},
    candidates::AbstractVector{<:Integer},
    X::AbstractMatrix,
    λ::Real;
    σ²a::Real = 1.0,
    w = nothing,
)
    nr = length(ref)
    size(X, 1) == nr || throw(DimensionMismatch("X needs one row per record in ref"))
    λ > 0 || throw(ArgumentError("λ must be positive"))
    V = Symmetric(Matrix{Float64}(K[ref, ref]) + λ * (w === nothing ? I : Diagonal(1 ./ w)))
    F = cholesky(V)
    Krc = Matrix{Float64}(K[ref, candidates])
    B = F \ Krc                                       # V⁻¹K_rc
    Xd = Matrix{Float64}(X)
    ViX = F \ Xd
    PB = B - ViX * (cholesky(Symmetric(Xd' * ViX)) \ (Xd' * B))
    Γ = Float64(σ²a) .* (Matrix{Float64}(K[candidates, candidates]) - Krc' * PB)
    (Γ .+ Γ') ./ 2
end

"""
    reliability(m::MME, k::Integer; zero = Int[], Kdiag = nothing) -> Matrix{Float64}

Reliability ``r^2 = 1 - \\operatorname{PEV}/\\sigma^2_u`` of the
levels of random effect `k`, per component (``\\sigma^2_u`` is the
diagonal of `G0` times the diagonal of ``\\mathbf{K}``). Without
`Kdiag`, ``\\operatorname{diag}(\\mathbf{K})`` is taken from the selected
inverse of a sparse ``\\mathbf{K}^{-1}`` (or a dense inverse of a dense
one); for ``\\mathbf{K} = A``, `Kdiag = nrm_diag(ped)` (i.e.
``1 + F_i``) is cheaper still. For non-inbred animals this is Henderson's
``1 - d_i\\alpha``.
"""
function reliability(m::MME, k::Integer; zero = Int[], Kdiag = nothing)
    r = m.effects[k]
    P = pev(m, k; zero = zero)
    kd = Kdiag === nothing ? _inverse_diagonal(r.Kinv) : Kdiag
    1 .- P ./ (kd * diag(r.G0)')
end

"""
    solutions(m::MME, x) -> (fixed = b, random = [u₁, u₂, …])

Split a solution vector into fixed effects and, per random effect, a
`q × m` matrix (levels in rows, components in columns).
"""
function solutions(m::MME, x::AbstractVector)
    fixed = x[m.blocks[1]]
    random = [
        reshape(x[m.blocks[k+1]], nlevels(r), ncomponents(r)) for
        (k, r) in enumerate(m.effects)
    ]
    (fixed = fixed, random = random)
end
