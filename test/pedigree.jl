@testset "Pedigree methods (RelationshipMatrices)" begin
    # inbred pedigree; parents are not the first rows
    ped = DataFrame(
        sire = [0, 0, 0, 1, 1, 4, 4, 0, 6, 6],
        dam = [0, 0, 0, 2, 3, 5, 5, 0, 7, 8],
    )
    rec = [4, 5, 6, 7, 9, 10]
    yy = [3.0, 2.5, 4.1, 3.3, 2.2, 3.9]
    α = 2.5
    F = nrm_diag(ped) .- 1
    @test maximum(F) > 0
    iod = AnimalModelIOD(yy, [ones(Int, 6)], rec, ped, α)
    m = mme(yy, ones(6, 1), [RandomEffect(incidence(rec, 10), ainv(ped), 1.0)]; σ²e = α)
    dd = randn(11)
    v = similar(dd)
    iod(v, dd)
    @test v ≈ m.lhs * dd
    @test mendelian_precision(ped) ≈ mendelian_precision(ped.sire, ped.dam; F)
    x = solve_mme(m)
    Wr, rinv, par, Api = ram_design(ped, rec, 1.0, α)
    @test par == [1, 2, 3, 4, 5, 6, 7, 8]
    mr = mme(yy, ones(6, 1), [RandomEffect(Wr, Api, 1.0)]; Rinv = rinv)
    b, (ap,) = solutions(mr, solve_mme(mr))
    @test b ≈ x[1:1]
    @test ram_backsolve(ped, rec, yy .- b[1], vec(ap), par, 1.0, α) ≈ x[2:end]
    # non-contiguous parents: animal 2 has no progeny
    ped2 = DataFrame(sire = [0, 0, 0, 1, 1, 4], dam = [0, 0, 0, 3, 3, 5])
    W2, _, par2, A2 = ram_design(ped2, [4, 5, 6], 1.0, 2.0)
    @test par2 == [1, 3, 4, 5]
    @test Matrix(A2) ≈ inv(nrm(ped2)[par2, par2])
    Mg = rand(Xoshiro(5), 0:2, 3, 4)
    r1 = base_allele_frequencies(Mg, ped2, [4, 5, 6])
    r2 = base_allele_frequencies(Mg, ainv(ped2), [4, 5, 6])
    @test r1[1] ≈ r2[1]
end
