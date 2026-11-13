@testset "Iteration on data and partitions" begin
    sexi = [1, 2, 2, 1, 1]
    iod = AnimalModelIOD(WWG, [sexi], CALF, SIRE, DAM, 2.0)
    m = mme(WWG, first(incidence(sexi)), [a]; σ²e = 40.0)
    d = randn(10)
    v = similar(d)
    iod(v, d)
    @test v ≈ m.lhs * d
    @test diag(iod) ≈ diag(m.lhs)
    @test rhs(iod) ≈ m.rhs
    x = solve_mme(m)
    @test first(gauss_seidel(iod; tol = 1e-26)) ≈ x atol = 1e-10
    @test first(pcg(iod, rhs(iod); M = diag(iod), tol = 1e-12)) ≈ x atol = 1e-8
    # n₁PA + n₂YD + n₃PC reproduces the solutions
    nrec = [i in CALF ? 1 : 0 for i = 1:8]
    yd = zeros(8)
    yd[CALF] = WWG - [x[1], x[2], x[2], x[1], x[1]]
    p = ebv_partition(SIRE, DAM, 2.0, nrec, yd, x[3:end])
    @test p.n1 .* p.PA + p.n2 .* p.YD + p.n3 .* replace(p.PC, NaN => 0) ≈ x[3:end]
    @test p.n1 + p.n2 + p.n3 ≈ ones(8)
    # multivariate partition with one trait equals the univariate one
    pm = ebv_partition_mt(
        SIRE,
        DAM,
        fill(20.0, 1, 1),
        [fill(n / 40, 1, 1) for n in nrec],
        reshape(yd, :, 1),
        reshape(x[3:end], :, 1),
    )
    @test [w[1] for w in pm.W1] ≈ p.n1
    @test [d[1] for d in pm.DYD][3] ≈ p.DYD[3]
end
