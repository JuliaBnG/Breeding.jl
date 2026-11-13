@testset "Reduced animal model equals the full model" begin
    Wr, rinv, par, _ = ram_design(SIRE, DAM, CALF, 20.0, 40.0)
    Ap = ainv_ref(SIRE[par], DAM[par])
    m = mme(WWG, X, [RandomEffect(Wr, sparse(Ap), 20.0)]; Rinv = rinv)
    b, (ap,) = solutions(m, solve_mme(m))
    full = solve_mme(mme(WWG, X, [a]; σ²e = 40.0))
    @test b ≈ full[1:2]
    abv = ram_backsolve(SIRE, DAM, CALF, WWG - X * b, vec(ap), par, 20.0, 40.0)
    @test abv ≈ full[3:end]
end
