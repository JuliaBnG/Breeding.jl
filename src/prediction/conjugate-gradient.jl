"""
    conjugate_gradient(A, b; max_iter=1000, tol=1e-6)

Solves the linear system Ax = b using the Conjugate Gradient iterative method.
This method is suitable for large, symmetric, positive-definite matrices.

# Arguments
- `A::AbstractMatrix`: The symmetric, positive-definite matrix.
- `b::AbstractVector`: The right-hand side vector.
- `max_iter::Int`: The maximum number of iterations.
- `tol::Float64`: The tolerance for the residual norm.

# Returns
- `x::AbstractVector`: The solution vector.
- `iterations::Int`: The number of iterations performed.
- `residual::Float64`: The final residual norm.
"""
function conjugate_gradient(
    A::Union{AbstractMatrix,Function},
    b::AbstractVector;
    max_iter::Int = 1000,
    tol::Real = 1e-6,
)
    # Problem size and element type
    n = length(b)
    T = eltype(b)

    # Preallocate working vectors to avoid allocations inside the loop
    x = zeros(T, n)
    r = copy(b)            # with x==0, r = b
    p = similar(r)
    copy!(p, r)
    Ap = similar(r)

    # Work with squared tolerance to avoid repeated sqrt calls
    r2 = LinearAlgebra.dot(r, r)
    tol2 = float(tol)^2

    if r2 <= tol2
        return x, 0, sqrt(r2)
    end

    iterations = 0
    for i = 1:max_iter
        iterations = i

        # Compute A * p into preallocated Ap without extra allocations.
        if A isa Function
            # Prefer in-place signature A!(p, out) if user provided one,
            # otherwise fall back to A(p) and copy the result.
            try
                A(p, Ap)
            catch _
                tmp = A(p)
                copy!(Ap, tmp)
            end
        else
            # mul! works efficiently for dense and sparse matrices
            mul!(Ap, A, p)
        end

        denom = LinearAlgebra.dot(p, Ap)
        if denom == 0
            # Breakdown: p is A-orthogonal to Ap; return current iterate
            return x, iterations, sqrt(r2)
        end

        alpha = r2 / denom

        # In-place updates to avoid temporaries
        @inbounds for j = 1:n
            x[j] += alpha * p[j]
            r[j] -= alpha * Ap[j]
        end

        r2_new = LinearAlgebra.dot(r, r)
        if r2_new <= tol2
            return x, iterations, sqrt(r2_new)
        end

        beta = r2_new / r2
        @inbounds for j = 1:n
            p[j] = r[j] + beta * p[j]
        end

        r2 = r2_new
    end

    return x, iterations, sqrt(r2)
end
