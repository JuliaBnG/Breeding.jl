@testset "Multi-component random effects" begin
    G0 = [20.0 -4; -4 10]
    W2 = incidence([2, 2, 2, 5, 6], 8)
    r = RandomEffect([Z, W2], sparse(Ai), G0)
    m = mme(WWG, X, [r]; σ²e = 40.0)
    W = [X Z W2]
    C = W'W + cat(zeros(2, 2), 40 .* kron(inv(G0), Ai); dims = (1, 2))
    @test Matrix(m.lhs) ≈ C
    x = solve_mme(m)
    @test x ≈ C \ (W'WWG)
    # maternal solutions are the second component, not the direct ones
    _, (u,) = solutions(m, x)
    @test u[:, 1] ≈ x[3:10] && u[:, 2] ≈ x[11:18]
    @test !(u[:, 2] ≈ u[:, 1])
    # warm start: x0 at the solution needs no iteration
    @test solve_mme(m; method = :pcg, x0 = x) == x
    @test solve_mme(m; method = :pcg, x0 = x, zero = [1]) ≈ solve_mme(m; zero = [1])
    @test_throws DimensionMismatch solve_mme(m; method = :pcg, x0 = x[1:3])
    # a dense system is read from one triangle only
    A = Matrix(m.lhs)
    A[1, 2:end] .= NaN
    @test first(pcg(A, m.rhs)) ≈ x
    @test first(pcg(permutedims(A), m.rhs; uplo = :U)) ≈ x
    @test first(pcg(permutedims(A), m.rhs; uplo = 'U')) ≈ x
end
