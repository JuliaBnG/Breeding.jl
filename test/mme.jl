@testset "MME, solvers and accuracies (Example 4.1)" begin
    m = mme(WWG, X, [a]; σ²e = 40.0)
    x = solve_mme(m)
    @test all(
        abs.(
            x -
            [4.358, 3.404, 0.098, -0.019, -0.041, -0.009, -0.186, 0.177, -0.249, 0.183],
        ) .< 6e-4,
    )
    # the R⁻¹ form gives the same solutions
    mr = mme(WWG, X, [a]; Rinv = fill(1 / 40, 5))
    @test solve_mme(mr) ≈ x
    # dense Henderson equations
    W = [X Z]
    C = W'W + cat(zeros(2, 2), 2Ai; dims = (1, 2))
    @test Matrix(m.lhs) ≈ C
    @test C \ (W'WWG) ≈ x
    for method in (:pcg, :gauss_seidel, :sor, :pinv)
        @test solve_mme(m; method, tol = 1e-26, maxiter = 10_000) ≈ x atol = 1e-10
    end
    @test solve_mme(m; method = :sor, ω = 1.1, tol = 1e-26) ≈ x atol = 1e-10
    @test first(jacobi(m.lhs, m.rhs; ω = [1; 1; fill(0.8, 8)], tol = 1e-26)) ≈ x atol =
        1e-10
    x0, it0, conv0 = pcg(sparse([2.0 0; 0 3]), zeros(2))
    @test x0 == zeros(2) && it0 == 0 && conv0 == 0
    xwarm, itwarm, convwarm = pcg(m.lhs, m.rhs; x0 = x)
    @test xwarm == x && itwarm == 0 && convwarm < 1e-12
    Cinv = lhs_inverse(m)
    @test vec(pev(m, 1)) ≈ diag(Cinv)[3:10] .* 40
    @test vec(reliability(m, 1)) ≈ 1 .- diag(Cinv)[3:10] .* 2
    b, (u,) = solutions(m, x)
    @test b ≈ x[1:2] && vec(u) ≈ x[3:10]
end
