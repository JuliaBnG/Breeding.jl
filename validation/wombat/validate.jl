# Validation of `Breeding.reml` against WOMBAT (Meyer, 2007).
#
#   julia -t 4 --project=<env with Breeding, RelationshipMatrices, DataFrames> \
#       validation/wombat/validate.jl [workdir]
#
# Two simulated data sets on one pedigree are written in WOMBAT's formats,
# analysed by `wombat` (which must be on the PATH) and by `Breeding.reml`
# from the same starting values, and the estimates and maximised
# restricted log-likelihoods are printed side by side.
#
# Simulation model (all stated, since this is the reference for the test):
#   - Pedigree: 6 discrete generations of 500 animals. Generation 1 are
#     unrelated, non-inbred founders. In later generations each animal has
#     a sire drawn at random from 25 males and a dam from 250 females of
#     the previous generation (sampling with replacement, so family sizes
#     are multinomial and inbreeding accumulates).
#   - Breeding values: founders N(0, G0); offspring ½(sire + dam) plus a
#     Mendelian sampling term with variance ½(1 − (F_s + F_d)/2) G0.
#   - Records on generations 2–6 (2500 animals); 20 contemporary groups
#     assigned at random, with effects N(0, 25) fitted as fixed.
#   - Analysis 1, one trait: σ²a = 30, σ²e = 70 (h² = 0.3).
#   - Analysis 2, two traits: G0 = [30 12; 12 20], R0 = [70 20; 20 60];
#     trait 2 is missing for a random 20% of the animals.
# Seeds are fixed below. One data set per analysis: this compares two
# programs on the same data, it is not a study of the estimator.
using Breeding, RelationshipMatrices, DataFrames, LinearAlgebra, SparseArrays, Random, Printf

const NGEN, NPER, NSIRE, NDAM, NHERD = 6, 500, 25, 250, 20

function simulate(; seed)
    rng = Xoshiro(seed)
    N = NGEN * NPER
    sire, dam = zeros(Int, N), zeros(Int, N)
    for i = (NPER+1):N
        lo = ((i - 1) ÷ NPER - 1) * NPER
        sire[i] = lo + rand(rng, 1:NSIRE)
        dam[i] = lo + NSIRE + rand(rng, 1:NDAM)
    end
    F = nrm_diag(sire, dam) .- 1
    G0 = [30.0 12.0; 12.0 20.0]
    R0 = [70.0 20.0; 20.0 60.0]
    Lg, Lr = cholesky(G0).L, cholesky(R0).L
    u = zeros(N, 2)
    for i = 1:N
        if sire[i] == 0
            u[i, :] = Lg * randn(rng, 2)
        else
            ms = sqrt(0.5 * (1 - (F[sire[i]] + F[dam[i]]) / 2))
            u[i, :] = 0.5 .* (u[sire[i], :] .+ u[dam[i], :]) .+ ms .* (Lg * randn(rng, 2))
        end
    end
    rec = collect((NPER+1):N)
    herd = rand(rng, 1:NHERD, length(rec))
    hfx = 5 .* randn(rng, NHERD, 2)
    y = [100.0 50.0] .+ hfx[herd, :] .+ u[rec, :] .+ (Lr * randn(rng, 2, length(rec)))'
    has2 = rand(rng, length(rec)) .> 0.2
    (; N, sire, dam, rec, herd, y, has2)
end

function write_wombat(dir, d)
    mkpath(joinpath(dir, "uni"))
    mkpath(joinpath(dir, "muv"))
    for sub in ("uni", "muv")
        open(joinpath(dir, sub, "pedigree.dat"), "w") do io
            for i = 1:d.N
                println(io, i, ' ', d.sire[i], ' ', d.dam[i])
            end
        end
    end
    open(joinpath(dir, "uni", "data.dat"), "w") do io
        for (k, a) in enumerate(d.rec)
            @printf(io, "%d %d %.6f\n", a, d.herd[k], d.y[k, 1])
        end
    end
    write(joinpath(dir, "uni", "wombat.par"), """
        COMMENT Breeding.jl validation: univariate animal model
        ANALYSIS UNI
        PEDS pedigree.dat
        DATA data.dat
          animal $(d.N)
          herd $NHERD
          y
        END
        MODEL
          FIX herd
          RAN animal NRM
          TRAIT y
        END MOD
        VAR animal 1
          20.0
        VAR residual 1
          80.0
        """)
    # multivariate: one line per record, trait number first, sorted by animal then trait
    open(joinpath(dir, "muv", "data.dat"), "w") do io
        for (k, a) in enumerate(d.rec)
            @printf(io, "1 %d %d %.6f\n", a, d.herd[k], d.y[k, 1])
            d.has2[k] && @printf(io, "2 %d %d %.6f\n", a, d.herd[k], d.y[k, 2])
        end
    end
    write(joinpath(dir, "muv", "wombat.par"), """
        COMMENT Breeding.jl validation: bivariate animal model, missing records
        ANALYSIS MUV 2
        PEDS pedigree.dat
        DATA data.dat GRP
          TRNOS 1 2
          traitno 2
          animal $(d.N)
          herd $NHERD
          y
        END
        MODEL
          FIX herd
          RAN animal NRM
          TRAIT y 1
          TRAIT y 2
        END MOD
        VAR animal 2
          20.0 5.0 15.0
        VAR residual 2
          80.0 10.0 50.0
        """)
end

function run_wombat(dir)
    cd(dir) do
        run(pipeline(`wombat wombat.par`; stdout = "wombat.stdout", stderr = "wombat.stderr"))
    end
    txt = read(joinpath(dir, "SumEstimates.out"), String)
    # "COVS A 1 2   estimate   s.e."
    cov(tag) = Dict((parse(Int, m[1]), parse(Int, m[2])) => (parse(Float64, m[3]), parse(Float64, m[4])) for
                    m in eachmatch(Regex("COVS $tag\\s+(\\d+)\\s+(\\d+)\\s+([-+.\\dEe]+)\\s+([-+.\\dEe]+)"), txt))
    m = match(r"Maximum log L\s*=?\s*([-+.\dEe]+)", read(joinpath(dir, "SumModel.out"), String) * txt)
    (A = cov("A"), Z = cov("Z"), logL = m === nothing ? NaN : parse(Float64, m[1]))
end

function fit_uni(d)
    X, _ = incidence(d.herd)
    Z = incidence(d.rec, d.N)
    Ai = ainv(d.sire, d.dam)
    r = reml(d.y[:, 1], X, [RandomEffect(Z, Ai, 20.0)]; σ²e = 80.0, tol = 1e-10)
    se = sqrt.(diag(r.AIinv))                  # order: residual, genetic
    (A = Dict((1, 1) => (r.G0[1][1, 1], se[2])), Z = Dict((1, 1) => (r.σ²e, se[1])), logL = r.logL,
        iterations = r.iterations)
end

function fit_muv(d)
    n = length(d.rec)
    keep2 = findall(d.has2)
    y = [d.y[:, 1]; d.y[keep2, 2]]
    trait = [fill(1, n); fill(2, length(keep2))]
    animal = [d.rec; d.rec[keep2]]
    herd = [d.herd; d.herd[keep2]]
    X, _ = incidence(collect(zip(trait, herd)))                 # herd within trait
    Za = incidence(animal, d.N)
    Zs = [Diagonal(Float64.(trait .== t)) * Za for t = 1:2]
    Ai = ainv(d.sire, d.dam)
    res = MultiTrait(trait, animal, [80.0 10.0; 10.0 50.0])     # one residual group per animal
    r = reml(y, X, [RandomEffect(Zs, Ai, [20.0 5.0; 5.0 15.0])]; residual = res, tol = 1e-10)
    G, R = r.G0[1], r.R0[1]
    se = sqrt.(diag(r.AIinv))                  # residual 11, 21, 22, then genetic 11, 21, 22
    keys3 = ((1, 1), (1, 2), (2, 2))
    (A = Dict(k => (G[k...], se[3+i]) for (i, k) in enumerate(keys3)),
        Z = Dict(k => (R[k...], se[i]) for (i, k) in enumerate(keys3)), logL = r.logL,
        iterations = r.iterations)
end

function report(name, w, b)
    println("\n== ", name)
    @printf("%-13s %12s %12s %10s %9s %9s\n", "parameter", "WOMBAT", "Breeding.jl", "rel. diff", "s.e. W", "s.e. B")
    for (tag, lab) in ((:A, "genetic"), (:Z, "residual")), k in sort(collect(keys(getfield(b, tag))))
        x, sx = getfield(w, tag)[k]
        v, sv = getfield(b, tag)[k]
        @printf("%-13s %12.4f %12.4f %10.1e %9.4f %9.4f\n", "$lab $(k[1]),$(k[2])", x, v, abs(v - x) / abs(x), sx, sv)
    end
    @printf("%-13s %12.3f %12.3f %10.1e   (absolute difference)\n", "log L", w.logL, b.logL, abs(b.logL - w.logL))
    println("Breeding.jl converged in ", b.iterations, " iterations")
end

dir = abspath(get(ARGS, 1, mktempdir()))
d = simulate(; seed = 20261113)
write_wombat(dir, d)
println("# ", strip(first(filter(l -> occursin("version", l), readlines(Cmd(`wombat --version`; dir = dir))))))
println("# Julia ", VERSION, ", ", Threads.nthreads(), " threads; work directory ", dir)
println("# ", d.N, " animals, ", length(d.rec), " with records, ", count(d.has2), " with trait 2")
report("Univariate animal model", run_wombat(joinpath(dir, "uni")), fit_uni(d))
report("Bivariate animal model, missing records", run_wombat(joinpath(dir, "muv")), fit_muv(d))
