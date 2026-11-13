@testset "Genomic helpers" begin
    rng = Xoshiro(11)
    M = rand(rng, 0:2, 12, 30)
    Zc, p, k = center_genotypes(M)
    @test all(abs.(sum(Zc; dims = 1)) .< 1e-10)
    G = Zc * Zc' / k + 0.01I
    y = randn(rng, 8)
    ref = 1:8
    λ = 3.0
    # SNP-BLUP and GBLUP agree; SNP effects are recovered from GBLUP
    msnp = mme(y, ones(8, 1), [RandomEffect(Zc[ref, :], I, 1 / k)]; σ²e = λ)
    mg = mme(
        y,
        ones(8, 1),
        [RandomEffect(incidence(ref, 12), inv(Zc * Zc' / k + 1e-9I), 1.0)];
        σ²e = λ,
    )
    xs, xg = solve_mme(msnp), solve_mme(mg)
    @test xs[1] ≈ xg[1] atol = 1e-6
    @test Zc * xs[2:end] ≈ xg[2:end] atol = 1e-5
    @test snp_from_gblup(Zc, inv(Zc * Zc' / k + 1e-9I), xg[2:end]; k) ≈ xs[2:end] atol =
        1e-5
    dgv, rel = gblup_index(G, ref, y .- xg[1], λ)
    @test length(dgv) == 12 && all(0 .≤ rel .≤ 1)
    H = rand(rng, 0:1, 10, 6)
    P, al = pseudo_snps(H, [1:3, 4:6])
    @test all(sum(P[:, 1:length(al[1])]; dims = 2) .== 2)
    p0, μ, dd = base_allele_frequencies(
        M[1:4, :],
        sparse(ainv_ref(zeros(Int, 6), zeros(Int, 6))),
        3:6,
    )
    @test μ ≈ vec(mean(M[1:4, :]; dims = 1)) atol = 0.05
end
