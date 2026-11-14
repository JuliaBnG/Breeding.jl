@testset "Selected inverse, PEV and reliability at scale" begin
    # a random pedigree of 400 animals with inbreeding, records on the last 300
    rng = Xoshiro(7)
    n = 400
    sire = zeros(Int, n)
    dam = zeros(Int, n)
    for i = 41:n
        sire[i] = rand(rng, 1:2:(i-1))
        dam[i] = rand(rng, 2:2:(i-1))
    end
    ped = DataFrame(sire = sire, dam = dam)
    Ainv = ainv(ped)
    rec = 101:n
    herd = rand(rng, 1:8, length(rec))
    Xh = sparse(1:length(rec), herd, 1.0, length(rec), 8)
    Zr = incidence(collect(rec), n)
    y = randn(rng, length(rec))
    m = mme(y, Xh, [RandomEffect(Zr, Ainv, 0.4)]; σ²e = 0.6)
    Cinv = lhs_inverse(m)

    # every stored element of the selected inverse equals the dense inverse,
    # and the pattern covers that of C
    S = selected_inverse(m)
    I_, J_, V_ = findnz(S)
    @test maximum(abs(V_[t] - Cinv[I_[t], J_[t]]) for t in eachindex(V_)) < 1e-10
    Ic, Jc, _ = findnz(m.lhs)
    @test all(S[Ic[t], Jc[t]] != 0 || abs(Cinv[Ic[t], Jc[t]]) < 1e-14 for t in eachindex(Ic))
    @test nnz(S) < length(Cinv)

    # PEV and reliability from the selected inverse against the dense inverse
    @test vec(pev(m, 1)) ≈ diag(Cinv)[9:end] .* 0.6
    @test diag(pev_candidates(m, 1, 350:400)) ≈ vec(pev(m, 1))[350:400]
    r2 = reliability(m, 1)
    @test vec(r2) ≈ 1 .- vec(pev(m, 1)) ./ (nrm_diag(ped) .* 0.4)
    @test reliability(m, 1; Kdiag = nrm_diag(ped)) ≈ r2
    @test all(-1e-10 .<= r2 .<= 1)
    # parent–offspring prediction error covariance is on the pattern
    @test S[8+n, 8+sire[n]] * 0.6 ≈ Cinv[8+n, 8+sire[n]] * 0.6

    # constrained equations: zeros there, the generalized inverse elsewhere
    Xo = [ones(length(rec)) Matrix(Xh)]                    # overparameterised
    mo = mme(y, Xo, [RandomEffect(Zr, Ainv, 0.4)]; σ²e = 0.6)
    So = selected_inverse(mo; zero = [1])
    Co = lhs_inverse(mo; zero = [1])
    @test all(iszero, So[1, :]) && all(iszero, So[:, 1])
    @test diag(So) ≈ diag(Co)
    @test vec(pev(mo, 1; zero = [1])) ≈ diag(Co)[10:end] .* 0.6
    @test vec(pev(mo, 1; zero = [1])) ≈ vec(pev(m, 1))
    @test_throws PosDefException selected_inverse(mo)
end
