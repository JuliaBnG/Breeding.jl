using Breeding
using DataFrames
using RelationshipMatrices
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using Test

# A⁻¹ from RelationshipMatrices for parent vectors
ainv_ref(s, d) = Matrix(ainv(DataFrame(sire = s, dam = d)))

# Mrode Example 4.1
const SIRE = [0, 0, 0, 1, 3, 1, 4, 3]
const DAM = [0, 0, 0, 0, 2, 2, 5, 6]
const CALF = [4, 5, 6, 7, 8]
const SEX = ["M", "F", "F", "M", "M"]
const WWG = [4.5, 2.9, 3.9, 3.5, 5.0]

@testset "Breeding: prediction and variance components" begin
    Ai = ainv_ref(SIRE, DAM)
    X, lv = incidence(SEX; levels = ["M", "F"])
    Z = incidence(CALF, 8)
    a = RandomEffect(Z, sparse(Ai), 20.0)

    @testset "Design helpers" begin
        @test lv == ["M", "F"]
        @test Matrix(X) == [1 0; 0 1; 0 1; 1 0; 1 0]
        @test incidence([2, 0, 1], 3) == sparse([1, 3], [2, 1], [1.0, 1.0], 3, 3)
        Φ = legendre([-1.0, 0.0, 1.0], 3; tmin = -1, tmax = 1)
        @test Φ[3, :] ≈ sqrt.([1, 3, 5] ./ 2)
        @test Φ[2, :] ≈ [sqrt(1 / 2), 0, -sqrt(5 / 2) / 2]
        F = linear_spline([4, 38, 106, 400], [4, 106, 208, 310])
        @test F[2, 1:2] ≈ [2 / 3, 1 / 3]
        @test F[3, :] ≈ [0, 1, 0, 0]
        @test F[4, 4] ≈ 310 / 400
        Zs = associate_incidence([1, 1, 1, 2, 2], [1, 2, 3, 4, 5], 5)
        @test Matrix(Zs) == [0 1 1 0 0; 1 0 1 0 0; 1 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0]
        Y = [1.0 missing; 2.0 3.0]
        y, obs = mt_stack(Y)
        @test y == [1.0, 2.0, 3.0]
        R0 = [2.0 0.5; 0.5 1.0]
        Ri = Matrix(mt_rinv(R0, obs))
        @test Ri[1, 1] ≈ 0.5
        @test Ri[[2, 3], [2, 3]] ≈ inv(R0)
        blocks = mt_blocks(obs, [ones(2, 1), ones(2, 1)])
        @test Matrix(reduce(hcat, blocks)) == [1 0; 1 0; 0 1]
    end

    @testset "MME, solvers and accuracies (Example 4.1)" begin
        m = mme(WWG, X, [a]; σ²e = 40.0)
        x = solve_mme(m)
        @test all(abs.(x - [4.358, 3.404, 0.098, -0.019, -0.041, -0.009, -0.186, 0.177,
                            -0.249, 0.183]) .< 6e-4)
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
        @test first(jacobi(m.lhs, m.rhs; ω = [1; 1; fill(0.8, 8)], tol = 1e-26)) ≈ x atol = 1e-10
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

    @testset "Multi-component random effects" begin
        G0 = [20.0 -4; -4 10]
        W2 = incidence([2, 2, 2, 5, 6], 8)
        r = RandomEffect([Z, W2], sparse(Ai), G0)
        m = mme(WWG, X, [r]; σ²e = 40.0)
        W = [X Z W2]
        C = W'W + cat(zeros(2, 2), 40 .* kron(inv(G0), Ai); dims = (1, 2))
        @test Matrix(m.lhs) ≈ C
        @test solve_mme(m) ≈ C \ (W'WWG)
    end

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
        nrec = [i in CALF ? 1 : 0 for i in 1:8]
        yd = zeros(8)
        yd[CALF] = WWG - [x[1], x[2], x[2], x[1], x[1]]
        p = ebv_partition(SIRE, DAM, 2.0, nrec, yd, x[3:end])
        @test p.n1 .* p.PA + p.n2 .* p.YD + p.n3 .* replace(p.PC, NaN => 0) ≈ x[3:end]
        @test p.n1 + p.n2 + p.n3 ≈ ones(8)
        # multivariate partition with one trait equals the univariate one
        pm = ebv_partition_mt(SIRE, DAM, fill(20.0, 1, 1), [fill(n / 40, 1, 1) for n in nrec],
                              reshape(yd, :, 1), reshape(x[3:end], :, 1))
        @test [w[1] for w in pm.W1] ≈ p.n1
        @test [d[1] for d in pm.DYD][3] ≈ p.DYD[3]
    end

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

    @testset "Pedigree methods (RelationshipMatrices)" begin
        # inbred pedigree; parents are not the first rows
        ped = DataFrame(sire = [0, 0, 0, 1, 1, 4, 4, 0, 6, 6], dam = [0, 0, 0, 2, 3, 5, 5, 0, 7, 8])
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
        mg = mme(y, ones(8, 1), [RandomEffect(incidence(ref, 12), inv(Zc * Zc' / k + 1e-9I),
                                              1.0)]; σ²e = λ)
        xs, xg = solve_mme(msnp), solve_mme(mg)
        @test xs[1] ≈ xg[1] atol = 1e-6
        @test Zc * xs[2:end] ≈ xg[2:end] atol = 1e-5
        @test snp_from_gblup(Zc, inv(Zc * Zc' / k + 1e-9I), xg[2:end]; k) ≈ xs[2:end] atol = 1e-5
        dgv, rel = gblup_index(G, ref, y .- xg[1], λ)
        @test length(dgv) == 12 && all(0 .≤ rel .≤ 1)
        H = rand(rng, 0:1, 10, 6)
        P, al = pseudo_snps(H, [1:3, 4:6])
        @test all(sum(P[:, 1:length(al[1])]; dims = 2) .== 2)
        p0, μ, dd = base_allele_frequencies(M[1:4, :], sparse(ainv_ref(zeros(Int, 6),
                                                                       zeros(Int, 6))), 3:6)
        @test μ ≈ vec(mean(M[1:4, :]; dims = 1)) atol = 0.05
    end

    @testset "Bayesian SNP models" begin
        rng = Xoshiro(3)
        Mg = Float64.(rand(rng, 0:2, 60, 15))
        g = zeros(15)
        g[[2, 7]] = [1.0, -1.0]
        y = 5 .+ Mg * g + 0.5randn(rng, 60)
        for method in (:A, :B, :C, :Cpi)
            r = bayes_snp(y, ones(60, 1), Mg; method, niter = 3000, burnin = 1000,
                          rng = Xoshiro(1))
            @test all(isfinite, r.g) && r.σ²e > 0
            @test r.g[2] > 0.5 && r.g[7] < -0.5
        end
    end

    @testset "Threshold models" begin
        # Mrode Example 15.1: first iteration from the book
        herd = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sex = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1]
        sire = [1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4]
        N = [1 0 0; 1 0 0; 1 0 0; 0 1 0; 1 0 1; 3 0 0; 1 1 0; 0 1 0; 1 0 0; 2 0 0; 1 0 0;
             0 0 1; 1 0 1; 1 0 0; 0 1 0; 0 0 1; 0 1 0; 1 0 0; 2 0 0; 2 0 0]
        Xt = [first(incidence(herd)) first(incidence(sex))]
        Zt = first(incidence(sire))
        r = threshold_model(N, Xt, Zt, ainv_ref([0, 0, 1, 3], [0, 0, 0, 0]), 19;
                            t0 = [0.468, 1.080], zero = [1, 3])
        @test r.history[:, 1] ≈ [0.441008, 1.044792, 0, 0.286869, 0, -0.358323, -0.041528,
                                 0.057853, 0.039850, -0.065178] atol = 1e-5
        @test round.([r.t; r.u], digits = 4) ≈ [0.4378, 1.0675, -0.0434, 0.0592, 0.0412,
                                                -0.0660]
        P = category_probabilities(r.t, [0.0, 1.0])
        @test sum(P; dims = 2) ≈ ones(2)
    end

    @testset "Covariance functions" begin
        G = [132.3 127.0 136.6; 127.0 172.8 200.8; 136.6 200.8 288.0]
        f = covariance_function(G, [90, 160, 240])
        @test f.Ĝ ≈ G
        @test f.M * f.T * f.M' ≈ G
        r = covariance_function(G, [90, 160, 240]; k = 2)
        @test size(r.C) == (2, 2) && r.df == 3
    end

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
        r1 = reml(y2, X, [RandomEffect(Z, sparse(Ai), r.G0[1][1])]; σ²e = r.σ²e, maxiter = 1)
        @test r1.logLs[1] ≈ Ld

        # Gibbs: with variances fixed, posterior means → BLUP
        blup = solve_mme(mme(WWG, X, [a]; σ²e = 40.0))
        means = reduce(hcat, [gibbs(WWG, X, [a]; σ²e = 40.0, estimate_variances = false,
                                    niter = 20_000, burnin = 500, rng = Xoshiro(s)).x
                              for s in 1:6])
        z = abs.(vec(mean(means; dims = 2)) - blup) ./
            (vec(std(means; dims = 2)) ./ sqrt(6))
        @test maximum(z) < 5
        g = gibbs(WWG, X, [a]; σ²e = 40.0, νe = 4, Se = 40.0, νu = [4],
                  Vu = [fill(80.0, 1, 1)], niter = 2000, burnin = 200, rng = Xoshiro(9))
        @test g.nsample == 1800 && all(>(0), g.samples)
        Y = [WWG WWG .+ randn(Xoshiro(2), 5)]
        gm = gibbs_mt(Y, X, Z, sparse(Ai); R0 = [40.0 10; 10 30], G0 = [20.0 5; 5 20],
                      νe = 10, Ve = [40.0 10; 10 30] * 7, νu = 10, Vu = [20.0 5; 5 20] * 7,
                      niter = 500, burnin = 100, rng = Xoshiro(4))
        @test size(gm.R) == (2, 2) && isposdef(gm.R) && isposdef(gm.G0[1])
    end
    @testset "REML: random regression, heterogeneous residuals" begin
        # simulated: 20 founders, 60 offspring with 6 records each on two
        # Legendre coefficients; genetic, PE, 3 residual classes, 12 groups
        rng = Xoshiro(2026)
        nf, no, nr = 20, 60, 6
        sire = [zeros(Int, nf); rand(rng, 1:2:nf, no)]
        dam = [zeros(Int, nf); rand(rng, 2:2:nf, no)]
        n = nf + no
        Ai = sparse(ainv(DataFrame(sire = sire, dam = dam)))
        A = inv(Matrix(Ai))
        animal = repeat(nf+1:n, inner = nr)
        t = repeat(1:nr, outer = no)
        Φ = legendre(t, 2)
        class = [x ≤ 2 ? 1 : x ≤ 4 ? 2 : 3 for x in t]
        grp = t .+ nr .* (animal .% 2)
        G0 = [4.0 1.0; 1.0 2.0]
        P0 = [3.0 0.5; 0.5 1.0]
        σ² = [2.0, 3.0, 5.0]
        Ag = cholesky(Symmetric(A)).L * randn(rng, n, 2) * cholesky(G0).U
        Pp = randn(rng, no, 2) * cholesky(P0).U
        y = [sum(Φ[r, c] * (Ag[animal[r], c] + Pp[animal[r]-nf, c]) for c in 1:2) + grp[r] +
             sqrt(σ²[class[r]]) * randn(rng) for r in eachindex(t)]
        X, _ = incidence(grp)
        Za = incidence(animal, n)
        Zp = incidence(animal .- nf, no)
        effs(G, P) = [RandomEffect([Diagonal(Φ[:, c]) * Za for c in 1:2], Ai, G),
                      RandomEffect([Diagonal(Φ[:, c]) * Zp for c in 1:2], I, P)]
        res = HetClasses(t, [1:2, 3:4, 5:6], [2.5, 2.5, 4.0])
        @test res.class == class
        @test_throws ArgumentError HetClasses(t, [1:2, 3:4], [1.0, 1.0])
        @test_throws ArgumentError HetClasses([1, 1, 3], [1.0, 1.0, 1.0])

        # log-likelihood equals its dense definition
        Gs, Ps, ss = [3.5 0.8; 0.8 1.5], [2.5 0.3; 0.3 1.2], [2.5, 2.5, 4.0]
        Zg = reduce(hcat, [Diagonal(Φ[:, c]) * Za for c in 1:2])
        Zq = reduce(hcat, [Diagonal(Φ[:, c]) * Zp for c in 1:2])
        V = Zg * kron(Gs, A) * Zg' + Zq * kron(Ps, Matrix(1.0I, no, no)) * Zq' +
            Diagonal(ss[class])
        Vi = inv(Symmetric(Matrix(V)))
        Xd = Matrix(X)                            # full column rank (one column per group)
        P = Vi - Vi * Xd * ((Xd' * Vi * Xd) \ (Xd' * Vi))
        Ld = -0.5 * (y' * P * y + logdet(Symmetric(Matrix(V))) + logdet(Symmetric(Xd' * Vi * Xd)))
        @test reml(y, X, effs(Gs, Ps); residual = res, maxiter = 0).logL ≈ Ld

        # AI (natural and Cholesky scale) and EM reach the same maximum
        fit(; kw...) = reml(y, X, effs(Gs, Ps); residual = HetClasses(class, ss), kw...)
        rn = fit(; tol = 1e-10)
        rc = fit(; param = :cholesky, tol = 1e-10)
        re = fit(; method = :em, tol = 1e-12, maxiter = 5000)
        θn = rn.history[:, end]
        @test rc.logL ≈ rn.logL atol = 1e-6
        @test rc.history[:, end] ≈ θn rtol = 1e-4
        @test re.logL ≈ rn.logL atol = 1e-3           # EM converges slowly
        @test re.history[:, end] ≈ θn rtol = 1e-2
        @test rn.residual == θn[1:3] && isnan(rn.σ²e)
        @test rn.G0[1] ≈ [θn[4] θn[5]; θn[5] θn[6]]
        @test size(rn.AIinv) == (9, 9) && isposdef(Symmetric(rn.AIinv))
        @test rn.iterations < 30 && rc.iterations < 30
        @test length(rc.steps) == rc.iterations && all(0 .≤ rc.steps .≤ 1)

        # stopping on log L change; ntol adds the expected gain of the next step
        rl = fit(; param = :cholesky, ltol = 0.5)
        rln = fit(; param = :cholesky, ltol = 0.5, ntol = 1e-6)
        @test rl.iterations ≤ rln.iterations
        @test rln.logL ≈ rc.logL atol = 1e-5

        # stationary point: central finite differences of logL vanish
        function L(θ)
            G = [θ[4] θ[5]; θ[5] θ[6]]
            Pm = [θ[7] θ[8]; θ[8] θ[9]]
            reml(y, X, effs(G, Pm); residual = HetClasses(class, θ[1:3]), maxiter = 0).logL
        end
        g = [(h = 1e-4 * max(1, abs(θn[i])); e = zeros(9); e[i] = h; (L(θn + e) - L(θn - e)) / 2h)
             for i in 1:9]
        @test maximum(abs.(g) .* max.(1, abs.(θn))) < 1e-3

        # sparse back-end (fixed pattern, CHOLMOD, selected inverse) = dense
        rs = fit(; tol = 1e-10, backend = :sparse)
        @test rs.logL ≈ rn.logL atol = 1e-8
        @test rs.history[:, end] ≈ θn rtol = 1e-6
        @test rs.AIinv ≈ rn.AIinv rtol = 1e-6
        @test rs.solutions ≈ rn.solutions rtol = 1e-6
        rs0 = reml(y, X, effs(Gs, Ps); residual = res, maxiter = 0, backend = :sparse)
        @test rs0.logL ≈ Ld
        # constrained fixed effects are dropped in the sparse back-end
        Xc = [ones(length(y)) X]                  # rank-deficient: intercept + groups
        rcd = reml(y, Xc, effs(Gs, Ps); residual = res, zero = [2], maxiter = 3, backend = :dense)
        rcs = reml(y, Xc, effs(Gs, Ps); residual = res, zero = [2], maxiter = 3, backend = :sparse)
        @test rcs.logLs ≈ rcd.logLs
        @test rcs.solutions ≈ rcd.solutions rtol = 1e-8

        # one class and Homogeneous agree with the σ²e form
        h1 = reml(y, X, effs(Gs, Ps); σ²e = 3.0, tol = 1e-10)
        h2 = reml(y, X, effs(Gs, Ps); residual = Homogeneous(3.0), tol = 1e-10)
        h3 = reml(y, X, effs(Gs, Ps); residual = HetClasses(ones(Int, length(y)), [3.0]),
                  tol = 1e-10)
        @test h1.history[:, end] ≈ h2.history[:, end]
        @test h3.logL ≈ h1.logL atol = 1e-6
        @test h1.σ²e ≈ h1.residual[1]
    end
end
