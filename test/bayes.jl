@testset "Bayesian SNP models" begin
    rng = Xoshiro(3)
    Mg = Float64.(rand(rng, 0:2, 60, 15))
    g = zeros(15)
    g[[2, 7]] = [1.0, -1.0]
    y = 5 .+ Mg * g + 0.5randn(rng, 60)
    for method in (:A, :B, :C, :Cpi)
        r = bayes_snp(
            y,
            ones(60, 1),
            Mg;
            method,
            niter = 3000,
            burnin = 1000,
            rng = Xoshiro(1),
        )
        @test all(isfinite, r.g) && r.σ²e > 0
        @test r.g[2] > 0.5 && r.g[7] < -0.5
    end
    # A SNP with z'z = 0 carries no information: its inclusion probability
    # is the prior 1 - π (BayesB, BayesC) or 1 - E(π | y) (BayesCπ), and its
    # effect is centred on zero
    Mz = copy(Mg)
    Mz[:, 15] .= 0
    for method in (:B, :C, :Cpi)
        r = bayes_snp(
            y,
            ones(60, 1),
            Mz;
            method,
            niter = 6000,
            burnin = 1000,
            rng = Xoshiro(2),
        )
        prior = method == :Cpi ? 1 - r.π : 0.7
        @test all(isfinite, r.g) && all(isfinite, r.inclusion)
        @test r.inclusion[15] ≈ prior atol = 0.05
        @test abs(r.g[15]) < 0.1
    end
end
