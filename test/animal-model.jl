@testset "animal_model" begin
    # Example 4.1 with labels, the pedigree in no particular order
    lab = ["a$i" for i = 1:8]
    ped = DataFrame(id = lab[[8, 3, 5, 1, 7, 2, 6, 4]])
    ped.sire = [s == 0 ? "" : lab[s] for s in SIRE[[8, 3, 5, 1, 7, 2, 6, 4]]]
    ped.dam = [d == 0 ? "" : lab[d] for d in DAM[[8, 3, 5, 1, 7, 2, 6, 4]]]
    dat = DataFrame(id = lab[CALF], sex = SEX, wwg = WWG,
                    dam = [d == 0 ? missing : lab[d] for d in DAM[CALF]])
    x41 = solve_mme(mme(WWG, X, [a]; σ²e = 40.0))     # X: M, F; animals 1:8

    @testset "single trait (Example 4.1)" begin
        m = animal_model(dat; traits = :wwg, fixed = :sex, random = [:id => (:A, 20.0)],
                         ped = ped, R0 = 40.0)
        @test m.residual isa Homogeneous && size(m.X) == (5, 2)
        e = blup(m; pev = true)
        u = e.random[1]
        @test names(u) == ["level", "id_wwg", "pev_id_wwg"]
        o = [findfirst(==(l), u.level) for l in lab]
        @test u.id_wwg[o] ≈ x41[3:10]
        f = e.fixed
        @test f.level[2:3] == ["F", "M"] && f.estimate[2] == 0.0
        @test f.estimate[1] ≈ x41[2] && f.estimate[3] ≈ x41[1] - x41[2]
        # PEV equal to those of the hand-built equations
        p41 = vec(pev(mme(WWG, X, [a]; σ²e = 40.0), 1))
        @test u.pev_id_wwg[o] ≈ p41
        # the same with the animals already numbered
        m2 = animal_model(DataFrame(id = CALF, sex = SEX, wwg = WWG); traits = [:wwg],
                          fixed = [:sex], random = [:id => (:A, 20.0)],
                          ped = DataFrame(id = 1:8, sire = SIRE, dam = DAM), R0 = 40.0)
        @test blup(m2).random[1].id_wwg ≈ x41[3:10]
        # a user K⁻¹ with labels is the same as :A
        m3 = animal_model(dat; traits = :wwg, fixed = :sex,
                          random = [:id => ((sparse(Ai), lab), 20.0)], R0 = 40.0)
        @test blup(m3).random[1].id_wwg ≈ x41[3:10]
    end

    @testset "direct and maternal effects" begin
        G0 = [20.0 -4; -4 10]
        m = animal_model(dat; traits = :wwg, fixed = :sex,
                         random = [(:id, :dam) => (:A, G0)], ped = ped, R0 = 40.0)
        W2 = incidence([0, 2, 2, 5, 6], 8)                # the dam of animal 4 is unknown
        x = solve_mme(mme(WWG, X, [RandomEffect([Z, W2], sparse(Ai), G0)]; σ²e = 40.0))
        u = blup(m).random[1]
        @test names(u) == ["level", "id_wwg", "dam_wwg"]
        o = [findfirst(==(l), u.level) for l in lab]
        @test u.id_wwg[o] ≈ x[3:10] && u.dam_wwg[o] ≈ x[11:18]
    end

    @testset "several traits, missing records" begin
        rng = Xoshiro(11)
        n = 60
        s = [0; 0; 0; 0; rand(rng, 1:4, n - 4)]
        d = [0; 0; 0; 0; 0; 0; rand(rng, 5:6, n - 6)]
        s[5:6] .= 0
        A1 = ainv(DataFrame(sire = s, dam = d))
        herd = rand(rng, ["h1", "h2", "h3"], n)
        age = randn(rng, n)
        y1 = randn(rng, n) .+ (herd .== "h2")
        y2 = Vector{Union{Missing,Float64}}(randn(rng, n))
        y2[rand(rng, 1:n, 15)] .= missing
        y1m = Vector{Union{Missing,Float64}}(y1)
        y1m[[3, 9]] .= missing
        dm = DataFrame(id = 1:n, herd = herd, age = age, y1 = y1m, y2 = y2)
        G0 = [1.0 0.3; 0.3 0.8]
        R0 = [2.0 0.5; 0.5 1.5]
        m = animal_model(dm; traits = [:y1, :y2], fixed = [:y1 => [:herd], :y2 => [:herd]],
                         covariates = Dict(:y2 => :age), random = [:id => (:A, G0)],
                         ped = DataFrame(id = 1:n, sire = s, dam = d), R0)
        # by hand: trait-major records, trait-specific means, herd and age effects
        Y = [y1m y2]
        yv, obs = mt_stack(Y)
        @test m.y == yv && length(m.y) == count(obs)
        H, _ = incidence(herd)
        Xs = [[ones(n) H[:, 2:3]], [ones(n) H[:, 2:3] age]]
        Xh = reduce(hcat, mt_blocks(obs, Xs))
        @test Matrix(m.X) == Matrix(Xh)
        Zh = mt_blocks(obs, [sparse(1.0I, n, n), sparse(1.0I, n, n)])
        eh = RandomEffect(Zh, A1, G0)
        Mh = mme(yv, Xh, [eh]; Rinv = mt_rinv(R0, obs))
        M = mme(m)
        @test M.lhs ≈ Mh.lhs && M.rhs ≈ Mh.rhs
        e = blup(m)
        @test e.solutions ≈ solve_mme(Mh)
        @test names(e.random[1]) == ["level", "id_y1", "id_y2"]
        @test nrow(e.fixed) == 2 * 4 + 1 && e.fixed.term[end] == :age
        # REML: the same likelihood as the hand-built model, and the estimates go back in
        r = reml(m; maxiter = 3)
        rh = reml(yv, Xh, [eh]; residual = MultiTrait(obs, R0), maxiter = 3)
        @test r.logL ≈ rh.logL && r.G0[1] ≈ rh.G0[1] && r.R0[1] ≈ rh.R0[1]
        m2 = animal_model(m, r)
        @test m2.effects[1].G0 ≈ r.G0[1] && m2.residual.R0[1] ≈ r.R0[1]
        @test mme(m2).lhs ≈ mme(yv, Xh, [RandomEffect(Zh, A1, r.G0[1])];
                                Rinv = mt_rinv(r.R0[1], obs)).lhs
        # an independent effect (permanent environment) over both traits
        dm.pe = dm.id
        mp = animal_model(dm; traits = [:y1, :y2], fixed = :herd,
                          random = [:id => (:A, G0), :pe => [0.5 0.1; 0.1 0.4]],
                          ped = DataFrame(id = 1:n, sire = s, dam = d), R0)
        q = length(unique(mp.row))                        # animals with a record
        @test mp.effects[2].Kinv == sparse(1.0I, q, q) && mp.levels[2] == sort(unique(mp.row))
        # one trait, REML through the model
        m1 = animal_model(dm; traits = :y1, fixed = :herd, random = [:id => (:A, 1.0)],
                          ped = DataFrame(id = 1:n, sire = s, dam = d), R0 = 2.0)
        r1 = reml(m1; maxiter = 3)
        o1 = .!ismissing.(y1m)
        r1h = reml(Float64.(y1m[o1]), sparse([ones(n) H[:, 2:3]][o1, :]),
                   [RandomEffect(incidence(findall(o1), n), A1, 1.0)]; σ²e = 2.0, maxiter = 3)
        @test r1.logL ≈ r1h.logL && r1.σ²e ≈ r1h.σ²e
        @test animal_model(m1, r1).residual.σ² == r1.σ²e
    end

    @testset "errors" begin
        P = DataFrame(id = 1:8, sire = SIRE, dam = DAM)
        d = DataFrame(id = CALF, sex = SEX, wwg = WWG, grp = SEX)
        kw = (; traits = :wwg, random = [:id => (:A, 20.0)], ped = P, R0 = 40.0)
        @test_throws ArgumentError animal_model(d; kw..., fixed = [:sex, :grp])  # confounded
        @test_throws ArgumentError animal_model(d; kw..., fixed = :nope)
        @test_throws ArgumentError animal_model(d; kw..., ped = nothing)
        @test_throws DimensionMismatch animal_model(d; kw..., random = [:id => (:A, [1.0 0; 0 1])])
        @test_throws DimensionMismatch animal_model(d; kw..., R0 = [1.0 0; 0 1])
        @test_throws ArgumentError animal_model(d; kw..., random = [:id => (:B, 1.0)])
        @test_throws ArgumentError animal_model(d; kw..., fixed = [:y9 => [:sex]])
        @test_throws ArgumentError animal_model(DataFrame(id = [9, 4], wwg = [1.0, 2.0]); kw...)
        @test_throws ArgumentError animal_model(DataFrame(id = [4], wwg = [missing]); kw...)
    end
end
