"""
    _conv(x, xold) -> Float64

Mrode's convergence criterion: the sum of squared changes in the
solutions divided by the sum of squared current solutions.
"""
function _conv(x, xold)
    num = den = 0.0
    @inbounds @simd for i in eachindex(x, xold)
        num += (x[i] - xold[i])^2
        den += x[i]^2
    end
    den == 0 ? (num == 0 ? 0.0 : Inf) : num / den
end

"""
    jacobi(A, b; x0 = zeros, ω = 1.0, second_order = false, tol = 1e-9,
           maxiter = 10_000, history = false) -> (x, iterations, conv[, hist])

Jacobi (total-step) iteration on ``\\mathbf{Ax} = \\mathbf{b}``, with
relaxation

```math
x^{(r+1)} = \\omega \\, M^{-1}(b - Ax^{(r)}) + x^{(r)},
```

where ``M = \\operatorname{diag}(A)``. `ω` is a scalar or a vector of
per-equation factors (Mrode uses 1 for fixed effects and 0.8 for
animal equations). With `second_order = true` the second-order Jacobi
``x^{(r+1)} = M^{-1}(b - Ax^{(r)}) + x^{(r)} + ω(x^{(r)} - x^{(r-1)})``
is used. Iteration stops when Mrode's convergence criterion (sum of squared
changes divided by sum of squared current solutions) is below `tol`. With
`history = true` the solutions of every round (columns, starting with
`x0`) are returned as a fourth value.
"""
function jacobi(A::AbstractMatrix, b::AbstractVector; x0 = zeros(length(b)), ω = 1.0,
                second_order = false, tol = 1e-9, maxiter = 10_000, history = false)
    n = length(b)
    x = Vector{Float64}(x0)
    xold = similar(x)
    xprev = copy(x)
    r = similar(x)
    d = Vector{Float64}(diag(A))
    any(iszero, d) && throw(ArgumentError("zero pivot on the diagonal"))
    hist = history ? [copy(x)] : Vector{Float64}[]
    conv = Inf
    it = 0
    while it < maxiter
        it += 1
        copyto!(xold, x)
        mul!(r, A, x)
        @inbounds for i in 1:n
            step = (b[i] - r[i]) / d[i]
            w = ω isa Number ? ω : ω[i]
            x[i] = second_order ? xold[i] + step + w * (xold[i] - xprev[i]) :
                   xold[i] + w * step
        end
        copyto!(xprev, xold)
        history && push!(hist, copy(x))
        conv = _conv(x, xold)
        conv < tol && break
    end
    history ? (x, it, conv, reduce(hcat, hist)) : (x, it, conv)
end

"""
    gauss_seidel(A, b; x0 = zeros, ω = 1.0, tol = 1e-9, maxiter = 10_000,
                 history = false) -> (x, iterations, conv[, hist])

Gauss–Seidel (ω = 1) or successive over-relaxation iteration on a
symmetric sparse system. Because ``\\mathbf{A}`` is symmetric, row `i`
is read from column `i` of the CSC storage, so one sweep costs one pass
over the non-zeros. Arguments as in [`jacobi`](@ref).
"""
function gauss_seidel(A::SparseMatrixCSC, b::AbstractVector; x0 = zeros(length(b)),
                      ω = 1.0, tol = 1e-9, maxiter = 10_000, history = false)
    n = length(b)
    x = Vector{Float64}(x0)
    xold = similar(x)
    rv, nz, cp = rowvals(A), nonzeros(A), A.colptr
    hist = history ? [copy(x)] : Vector{Float64}[]
    conv = Inf
    it = 0
    while it < maxiter
        it += 1
        copyto!(xold, x)
        @inbounds for j in 1:n
            s = b[j]
            dj = 0.0
            for p in cp[j]:cp[j+1]-1
                i = rv[p]
                if i == j
                    dj = nz[p]
                else
                    s -= nz[p] * x[i]
                end
            end
            dj == 0 && throw(ArgumentError("zero pivot in equation $j"))
            x[j] += ω * (s / dj - x[j])
        end
        history && push!(hist, copy(x))
        conv = _conv(x, xold)
        conv < tol && break
    end
    history ? (x, it, conv, reduce(hcat, hist)) : (x, it, conv)
end

gauss_seidel(A::AbstractMatrix, b::AbstractVector; kw...) = gauss_seidel(sparse(A), b; kw...)

"""
    pcg(A, b; M = diag(A), x0 = zeros, tol = 1e-12, maxiter = length(b) * 10,
        criterion = :residual, history = false) -> (x, iterations, conv[, hist])

Preconditioned conjugate gradient for a symmetric positive (semi-)
definite system (Lidauer et al., 1999; Mrode Section 19.6). `A` is a
matrix or a function `A(v, d)` that overwrites `v` with ``\\mathbf{A}d``
– e.g. an iteration-on-data operator such as
[`AnimalModelIOD`](@ref) – in which case `M` (the diagonal
preconditioner) must be given. Convergence is ``\\|b - Ax\\| / \\|b\\|``
(`criterion = :residual`) or Mrode's change criterion
(`criterion = :change`).
"""
function pcg(A, b::AbstractVector; M = nothing, x0 = zeros(length(b)), tol = 1e-12,
             maxiter = 10length(b), criterion::Symbol = :residual, history = false)
    n = length(b)
    Minv = 1 ./ Vector{Float64}(M === nothing ? diag(A) : M)
    x = Vector{Float64}(x0)
    xold = copy(x)
    v = similar(x)
    _amul!(v, A, x)
    e = b .- v                 # residual
    z = Minv .* e
    d = copy(z)
    ez = dot(e, z)
    nb = norm(b)
    hist = history ? [copy(x)] : Vector{Float64}[]
    conv = nb == 0 ? norm(e) : norm(e) / nb
    it = 0
    conv < tol && return history ? (x, it, conv, reduce(hcat, hist)) : (x, it, conv)
    while it < maxiter
        it += 1
        _amul!(v, A, d)
        w = ez / dot(d, v)
        copyto!(xold, x)
        @. x += w * d
        @. e -= w * v
        history && push!(hist, copy(x))
        conv = criterion == :change ? _conv(x, xold) : (nb == 0 ? norm(e) : norm(e) / nb)
        conv < tol && break
        @. z = Minv * e
        ez_new = dot(e, z)
        β = ez_new / ez
        ez = ez_new
        @. d = z + β * d
    end
    history ? (x, it, conv, reduce(hcat, hist)) : (x, it, conv)
end

_amul!(v, A::AbstractMatrix, d) = mul!(v, A, d)
_amul!(v, A, d) = A(v, d)
