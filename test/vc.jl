@testset "Variance components" begin
    # method 3 and canonical sums of squares (Mrode Section 17.3)
    y = [2.9, 4.0, 3.5, 3.5]
    Zs = Matrix(first(incidence([2, 1, 3, 2])))
    h = henderson3(y, ones(4, 1), Zs)
    @test [h.F, h.S, h.R] ≈ [48.3025, 0.4275, 0.18]
    c = sire_contrasts(y, ones(4, 1), Zs, Matrix(1.0I, 3, 3))
    @test c.w ≈ [1.5, 1.0]
    θ, V = fit_mean_squares([c.ss; c.R / c.dfR], [1, 1, c.dfR], [ones(3) [c.w; 0]])
    @test round.(θ, digits = 3) ≈ [0.163, 0.047]

    # REML (Mrode Section 17.7): AI iterates, logL, EM agrees
    y2 = [2.6, 0.1, 1.0, 3.0, 1.0]
    a2 = RandomEffect(Z, sparse(Ai), 0.2)
    r = reml(y2, X, [a2]; σ²e = 0.4)
    @test round.(r.history[:, 2], digits = 4) ≈ [0.4838, 0.3695]
    @test round(r.logLs[1], digits = 4) ≈ -2.3852
    @test [r.σ²e, r.G0[1][1]] ≈ [0.4835, 0.5514] atol = 1e-4
    em = reml(y2, X, [a2]; σ²e = 0.4, method = :em, maxiter = 5000, tol = 1e-12)
    @test [em.σ²e, em.G0[1][1]] ≈ [r.σ²e, r.G0[1][1]] atol = 1e-3
    # the REML log-likelihood equals its dense definition
    A = inv(Ai)
    V = Z * A * Z' .* r.G0[1][1] + r.σ²e * I
    Vi = inv(V)
    P = Vi - Vi * X * ((X' * Vi * X) \ (X' * Vi))
    Ld = -0.5 * (y2' * P * y2 + logdet(V) + logdet(X' * Vi * X))
    r1 =
        reml(y2, X, [RandomEffect(Z, sparse(Ai), r.G0[1][1])]; σ²e = r.σ²e, maxiter = 1)
    @test r1.logLs[1] ≈ Ld

    # Gibbs: with variances fixed, posterior means → BLUP
    blup = solve_mme(mme(WWG, X, [a]; σ²e = 40.0))
    means = reduce(
        hcat,
        [
            gibbs(
                WWG,
                X,
                [a];
                σ²e = 40.0,
                estimate_variances = false,
                niter = 20_000,
                burnin = 500,
                rng = Xoshiro(s),
            ).x for s = 1:6
        ],
    )
    z =
        abs.(vec(mean(means; dims = 2)) - blup) ./
        (vec(std(means; dims = 2)) ./ sqrt(6))
    @test maximum(z) < 5
    g = gibbs(
        WWG,
        X,
        [a];
        σ²e = 40.0,
        νe = 4,
        Se = 40.0,
        νu = [4],
        Vu = [fill(80.0, 1, 1)],
        niter = 2000,
        burnin = 200,
        rng = Xoshiro(9),
    )
    @test g.nsample == 1800 && all(>(0), g.samples)
    Y = [WWG WWG .+ randn(Xoshiro(2), 5)]
    gm = gibbs_mt(
        Y,
        X,
        Z,
        sparse(Ai);
        R0 = [40.0 10; 10 30],
        G0 = [20.0 5; 5 20],
        νe = 10,
        Ve = [40.0 10; 10 30] * 7,
        νu = 10,
        Vu = [20.0 5; 5 20] * 7,
        niter = 500,
        burnin = 100,
        rng = Xoshiro(4),
    )
    @test size(gm.R) == (2, 2) && isposdef(gm.R) && isposdef(gm.G0[1])
end
