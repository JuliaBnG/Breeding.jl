@testset "mtpblup" begin
    ped = DataFrame(sire = SIRE, dam = DAM)
    ped.ft_a = Union{Missing,Float64}[missing, missing, missing, WWG...]
    ped.ebv_a = zeros(8)
    Ais = sparse(Ai)
    # one trait: h² = 1/3, σ²a = 20, so σ²e = 40 as in Example 4.1
    mtpblup(ped, [Trait("a"; h² = 1 / 3)], Ais, fill(20.0, 1, 1))
    m1 = mme(WWG, ones(5, 1), [RandomEffect(Z, Ais, 20.0)]; σ²e = 40.0)
    @test ped.ebv_a ≈ solve_mme(m1)[2:end] atol = 1e-8
    # two traits: the equations of mme with R⁻¹ = diag(1/σ²e), σ²e = G_tt(1-h²)/h²
    # (a non-symmetric G⁻¹Λ penalty used to be fed to CG)
    ped.ft_b =
        Union{Missing,Float64}[missing, missing, missing, 1.0, 2.5, missing, 0.5, 3.0]
    ped.ebv_b = zeros(8)
    trts = [Trait("a"; h² = 1 / 3), Trait("b"; h² = 0.25)]
    G = [20.0 5; 5 10]
    p2 = copy(ped)
    mtpblup(ped, trts, Ais, G)
    mtpblup(p2, trts, Ais, G; method = :direct)
    @test ped.ebv_a ≈ p2.ebv_a atol = 1e-8
    @test ped.ebv_b ≈ p2.ebv_b atol = 1e-8
    ob = [4, 5, 7, 8] # both traits observed
    Zo = incidence(ob, 8)
    Zs2 = [[Zo; spzeros(4, 8)], [spzeros(4, 8); Zo]]
    Xo = kron(I(2), ones(4))
    σ²e = [40.0, 30.0]
    mm = mme(
        [ped.ft_a[ob]; ped.ft_b[ob]],
        Xo,
        [RandomEffect(Zs2, Ais, G)];
        Rinv = repeat(1 ./ σ²e; inner = 4),
    )
    @test [p2.ebv_a; p2.ebv_b] ≈ solve_mme(mm)[3:end]
    # auto-creation of EBV columns when absent, and validation on empty records
    p_no_ebv = select(ped, Not([:ebv_a, :ebv_b]))
    mtpblup(p_no_ebv, trts, Ais, G)
    @test p_no_ebv.ebv_a ≈ ped.ebv_a atol = 1e-8
    @test p_no_ebv.ebv_b ≈ ped.ebv_b atol = 1e-8
    p_empty = copy(ped)
    p_empty.ft_a .= missing
    @test_throws ArgumentError mtpblup(p_empty, trts, Ais, G)
    # a missing EBV column does not discard the other trait's warm start
    p_half = select(ped, Not(:ebv_b))
    mtpblup(p_half, trts, Ais, G)
    @test p_half.ebv_a ≈ ped.ebv_a atol = 1e-8
    @test p_half.ebv_b ≈ ped.ebv_b atol = 1e-8
end
