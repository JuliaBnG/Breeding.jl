@testset "OCS and Robust OCS (Meuwissen 1997; Fogg et al. 2026)" begin
    # Example 4.1 MME solutions and prediction error covariance Ω
    m_ocs = mme(WWG, X, [a]; σ²e = 40.0)
    sol_ocs = solve_mme(m_ocs)
    _, (u_mat,) = solutions(m_ocs, sol_ocs)
    u_ebv = vec(u_mat)
    Cinv_ocs = lhs_inverse(m_ocs)
    Ω_pec = Matrix(Cinv_ocs[3:10, 3:10]) .* 40.0
    Σ_rel = nrm(DataFrame(sire = SIRE, dam = DAM))
    sex_cand = [1, 0, 1, 1, 0, 0, 1, 1]

    # Candidate PEV: mixed-model form, V form and the full inverse agree
    Γ_mme = pev_candidates(m_ocs, 1, 1:8)
    @test Γ_mme ≈ Ω_pec
    @test diag(Γ_mme) ≈ vec(pev(m_ocs, 1))
    @test issymmetric(Γ_mme)
    Γ_v = pev_candidates(Σ_rel, CALF, 1:8, Matrix(X), 40.0 / 20.0; σ²a = 20.0)
    @test Γ_v ≈ Ω_pec
    sub = [8, 2, 5]
    @test pev_candidates(m_ocs, 1, sub) ≈ Ω_pec[sub, sub]
    @test pev_candidates(Σ_rel, CALF, sub, Matrix(X), 2.0; σ²a = 20.0) ≈ Ω_pec[sub, sub]
    @test_throws ArgumentError pev_candidates(m_ocs, 1, [9])
    @test_throws ArgumentError pev_candidates(m_ocs, 1, [1]; component = 2)
    @test_throws DimensionMismatch pev_candidates(Σ_rel, CALF, sub, Matrix(X)[1:3, :], 2.0)
    r_piped = rocs(u_ebv, Σ_rel; Ω = Γ_v, sex = sex_cand, λ = 2.0, κ = 1.0, ub = 0.25)
    r_full = rocs(u_ebv, Σ_rel; Ω = Ω_pec, sex = sex_cand, λ = 2.0, κ = 1.0, ub = 0.25)
    @test r_piped.w ≈ r_full.w atol = 1e-8

    # 1. TM1997 (DataFrame and vector methods agree and satisfy constraints)
    ped_ocs = DataFrame(tmi = u_ebv, sex = sex_cand)
    c_df = TM1997(ped_ocs, Σ_rel, 0.35)
    c_vec = TM1997(u_ebv, sex_cand, Σ_rel, 0.35)
    @test c_df ≈ c_vec
    @test all(>=(-1e-12), c_df)
    @test sum(c_df[sex_cand .== 1]) ≈ 0.5 atol = 1e-10
    @test sum(c_df[sex_cand .== 0]) ≈ 0.5 atol = 1e-10
    @test c_df' * Σ_rel * c_df ≤ 0.35 + 1e-8

    # OC_fixsex2: male optimization given fixed female contributions
    males_idx = findall(==(1), sex_cand)
    females_idx = findall(==(0), sex_cand)
    ord_idx = [males_idx; females_idx]
    c_m = OC_fixsex2(
        u_ebv[males_idx],
        Σ_rel[ord_idx, ord_idx],
        0.40,
        ones(length(females_idx)),
    )
    @test length(c_m) == length(males_idx)
    @test sum(c_m) ≈ 1.0 atol = 1e-8
    @test all(>=(-1e-12), c_m)

    # 2. Fogg et al. (2026) 4 regimes: TS (0,0), RTS (0,κ), OCS (λ,0), ROCS (λ,κ)
    r_ts = rocs(u_ebv, Σ_rel; Ω = Ω_pec, sex = sex_cand, λ = 0.0, κ = 0.0, ub = 0.25)
    r_rts = rocs(u_ebv, Σ_rel; Ω = Ω_pec, sex = sex_cand, λ = 0.0, κ = 1.0, ub = 0.25)
    r_ocs = rocs(u_ebv, Σ_rel; Ω = Ω_pec, sex = sex_cand, λ = 2.0, κ = 0.0, ub = 0.25)
    r_rocs = rocs(u_ebv, Σ_rel; Ω = Ω_pec, sex = sex_cand, λ = 2.0, κ = 1.0, ub = 0.25)

    for r in (r_ts, r_rts, r_ocs, r_rocs)
        @test r.converged
        @test sum(r.w[sex_cand .== 1]) ≈ 0.5 atol = 1e-8
        @test sum(r.w[sex_cand .== 0]) ≈ 0.5 atol = 1e-8
        @test all(-1e-12 .≤ r.w .≤ 0.25 + 1e-12)
        @test r.merit ≈ dot(r.w, u_ebv)
        @test r.coancestry ≈ dot(r.w, Σ_rel * r.w)
        @test r.std_error ≈ sqrt(dot(r.w, Ω_pec * r.w))
        @test r.objective ≈ r.merit - r.κ * r.std_error - 0.5 * r.λ * r.coancestry
    end

    # Qualitative ordering from Fogg et al. (2026) Table 2:
    # - TS maximizes expected genetic merit
    # - RTS reduces std_error relative to TS
    # - OCS reduces coancestry relative to TS
    # - ROCS reduces std_error relative to OCS and coancestry relative to RTS
    @test r_ts.merit ≥ r_rts.merit - 1e-9
    @test r_ts.merit ≥ r_ocs.merit - 1e-9
    @test r_rts.std_error < r_ts.std_error
    @test r_ocs.coancestry < r_ts.coancestry
    @test r_rocs.std_error < r_ocs.std_error
    @test r_rocs.coancestry < r_rts.coancestry

    # Tuple destructuring and vector PEV input
    w_d, obj_d, merit_d, coanc_d, se_d =
        robust_ocs(u_ebv, Σ_rel; Ω = diag(Ω_pec), sex = sex_cand, λ = 1.0, κ = 0.5)
    @test length(w_d) == 8 && isfinite(obj_d) && merit_d > 0 && coanc_d > 0 && se_d > 0
end
