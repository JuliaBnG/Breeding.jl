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

    # MultiTrait with one trait and one record per group = HetClasses
    m1 = reml(y, X, effs(Gs, Ps); residual = MultiTrait(ones(Int, length(y)), 1:length(y),
              [fill(s, 1, 1) for s in ss]; class = class), maxiter = 2)
    h4 = reml(y, X, effs(Gs, Ps); residual = HetClasses(class, ss), maxiter = 2)
    @test m1.logLs ≈ h4.logLs && m1.residual ≈ h4.residual
end
@testset "REML: multi-trait random regression, missing records" begin
    # two traits recorded at times 1:6 on 60 offspring, about a quarter
    # of the records missing; genetic RR (Legendre(2) per trait, 4 × 4
    # G0), PE intercept per trait, residual 2 × 2 per animal and time
    # with two classes of time
    rng = Xoshiro(11)
    nf, no, nt = 20, 60, 6
    sire = [zeros(Int, nf); rand(rng, 1:2:nf, no)]
    dam = [zeros(Int, nf); rand(rng, 2:2:nf, no)]
    n = nf + no
    Ai = sparse(ainv(DataFrame(sire = sire, dam = dam)))
    A = inv(Matrix(Ai))
    G0 = [4.0 1.0 2.0 0.5; 1.0 2.0 0.5 0.8; 2.0 0.5 5.0 1.0; 0.5 0.8 1.0 2.5]
    P0 = [3.0 1.0; 1.0 2.0]
    R0 = [[4.0 1.5; 1.5 3.0], [6.0 2.0; 2.0 5.0]]
    Ag = cholesky(Symmetric(A)).L * randn(rng, n, 4) * cholesky(G0).U
    Pp = randn(rng, no, 2) * cholesky(P0).U
    animal, time, trait, y = Int[], Int[], Int[], Float64[]
    Φt = legendre(1:nt, 2)
    for a in nf+1:n, s in 1:nt
        c = s ≤ 3 ? 1 : 2
        e = cholesky(R0[c]).L * randn(rng, 2)
        for tr in 1:2
            (tr == 2 ? rand(rng) < 0.35 : rand(rng) < 0.1) && continue
            push!(animal, a); push!(time, s); push!(trait, tr)
            u = sum(Φt[s, k] * Ag[a, 2(tr-1)+k] for k in 1:2)
            push!(y, 10tr + 0.5s + u + Pp[a-nf, tr] + e[tr])
        end
    end
    class = [s ≤ 3 ? 1 : 2 for s in time]
    Φ = legendre(time, 2)
    Za = incidence(animal, n)
    Zp = incidence(animal .- nf, no)
    X, _ = incidence(2 .* time .+ trait)     # trait × time means
    Zg = [Diagonal(Φ[:, k] .* (trait .== tr)) * Za for tr in 1:2 for k in 1:2]
    Zq = [Diagonal(Float64.(trait .== tr)) * Zp for tr in 1:2]
    effs(G, P) = [RandomEffect(Zg, Ai, G), RandomEffect(Zq, I, P)]
    grp = collect(zip(animal, time))
    mt(Rs) = MultiTrait(trait, grp, Rs; class = class)
    @test_throws ArgumentError MultiTrait([1, 1], [1, 1], [1.0 0; 0 1])      # trait twice
    @test_throws ArgumentError reml(y, X, effs(G0, P0);                     # never together
                                    residual = MultiTrait(trait, 1:length(y), [4.0 1; 1 3]))

    # log-likelihood equals its dense definition
    Gs = G0 + 0.5I
    Ps = [2.5 0.5; 0.5 1.5]
    Rs = [[3.0 1.0; 1.0 2.5], [5.0 1.0; 1.0 4.0]]
    Zgm = reduce(hcat, Zg)
    Zqm = reduce(hcat, Zq)
    Rd = zeros(length(y), length(y))
    for i in eachindex(y), j in eachindex(y)
        grp[i] == grp[j] && (Rd[i, j] = Rs[class[i]][trait[i], trait[j]])
    end
    V = Zgm * kron(Gs, A) * Zgm' + Zqm * kron(Ps, Matrix(1.0I, no, no)) * Zqm' + Rd
    Vi = inv(Symmetric(Matrix(V)))
    Xd = Matrix(X)
    P = Vi - Vi * Xd * ((Xd' * Vi * Xd) \ (Xd' * Vi))
    Ld = -0.5 * (y' * P * y + logdet(Symmetric(Matrix(V))) + logdet(Symmetric(Xd' * Vi * Xd)))
    @test reml(y, X, effs(Gs, Ps); residual = mt(Rs), maxiter = 0).logL ≈ Ld
    @test reml(y, X, effs(Gs, Ps); residual = mt(Rs), maxiter = 0, backend = :sparse).logL ≈ Ld

    # AI (natural and Cholesky scale) and EM reach the same maximum
    fit(; kw...) = reml(y, X, effs(Gs, Ps); residual = mt(Rs), kw...)
    rn = fit(; tol = 1e-10)
    rc = fit(; param = :cholesky, tol = 1e-10)
    re = fit(; method = :em, tol = 1e-12, maxiter = 5000)
    θn = rn.history[:, end]
    np = length(θn)
    @test np == 6 + 10 + 3
    @test rc.logL ≈ rn.logL atol = 1e-6
    @test rc.history[:, end] ≈ θn rtol = 1e-4
    @test re.logL ≈ rn.logL atol = 1e-3
    @test all(diff(re.logLs) .> -1e-8)                 # EM never decreases log L
    @test rn.R0[2] ≈ [θn[4] θn[5]; θn[5] θn[6]] && rn.residual == θn[1:6]
    @test isposdef(Symmetric(rn.AIinv)) && rn.iterations < 30

    # stationary point: central finite differences of logL vanish
    function L(θ)
        m = Breeding._unpack(θ, [2, 2, 4, 2])
        reml(y, X, effs(m[3], m[4]); residual = mt(m[1:2]), maxiter = 0).logL
    end
    g = [(h = 1e-4 * max(1, abs(θn[i])); e = zeros(np); e[i] = h; (L(θn + e) - L(θn - e)) / 2h)
         for i in 1:np]
    @test maximum(abs.(g) .* max.(1, abs.(θn))) < 1e-3

    # sparse back-end = dense
    rs = fit(; tol = 1e-10, backend = :sparse)
    @test rs.logL ≈ rn.logL atol = 1e-8
    @test rs.history[:, end] ≈ θn rtol = 1e-6
    @test rs.AIinv ≈ rn.AIinv rtol = 1e-6

    # constraints: a residual covariance fixed, two variances equal;
    # the constrained maximum is stationary in the free directions
    rf = fit(; fix = [5], equal = [[1, 4]], tol = 1e-10)
    θf = rf.history[:, end]
    @test all(rf.history[5, :] .== Rs[2][1, 2]) && θf[1] == θf[4]
    @test rf.logL < rn.logL
    dirs = [Matrix(1.0I, np, np)[:, setdiff(1:np, [1, 4, 5])] [1; zeros(2); 1; zeros(np - 4)]]
    gd = [(h = 1e-4 * max(1, maximum(abs, d .* θf)); (L(θf + h * d) - L(θf - h * d)) / 2h)
          for d in eachcol(dirs)]
    @test maximum(abs, gd) < 1e-3
    @test size(rf.AIinv) == (np, np) && all(rf.AIinv[5, :] .== 0)
    @test_throws ArgumentError fit(; fix = [5], param = :cholesky)
    @test_throws ArgumentError fit(; fix = [5], equal = [[5, 6]])

    # wide layout: rows of an n × t table are the residual groups
    Yw = [missing 1.0; 2.0 3.0; 4.0 missing; 5.0 6.5; 7.0 8.0; 2.5 missing]
    yw, ow = mt_stack(Yw)
    mw = MultiTrait(ow, [2.0 0.5; 0.5 1.0])
    @test mw.trait == [1, 1, 1, 1, 1, 2, 2, 2, 2] && mw.group == [1, 2, 3, 4, 5, 6, 1, 3, 4]
    Xw = reduce(hcat, mt_blocks(ow, [ones(6, 1), ones(6, 1)]))
    Vw = Matrix(inv(Matrix(mt_rinv([2.0 0.5; 0.5 1.0], ow))))
    Viw = inv(Symmetric(Vw))
    Pw = Viw - Viw * Xw * ((Xw' * Viw * Xw) \ (Xw' * Viw))
    Lw = -0.5 * (yw' * Pw * yw + logdet(Symmetric(Vw)) + logdet(Symmetric(Xw' * Viw * Xw)))
    @test reml(yw, Xw, RandomEffect[]; residual = mw, maxiter = 0).logL ≈ Lw
    @test nested([1.0 2; 3 4; 5 6], [2, 1, 2]) == sparse([0 0 1 2; 3 4 0 0; 0 0 5 6])
end
