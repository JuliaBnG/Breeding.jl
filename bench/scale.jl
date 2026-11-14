# Scale benchmark for the sparse paths of Breeding.jl.
#
#   julia -t 8 --project=bench bench/scale.jl 10000 100000
#
# Model: single-trait animal model, y = Xb + Za + e, on a simulated pedigree.
# Pedigree: discrete generations of `N ÷ 10` animals; each animal has a sire
# drawn from 5% and a dam from 50% of the previous generation (random
# mating among the selected, so inbreeding accumulates). Every animal after
# the first generation has one record; 50 contemporary groups; h² = 0.3.
# Not a validation: it measures time and memory, and checks PEV against
# solved columns for a sample of animals.
using Breeding, RelationshipMatrices, DataFrames, SparseArrays, LinearAlgebra, Random, Printf

function pedigree(N; rng)
    g = N ÷ 10
    sire, dam = zeros(Int, N), zeros(Int, N)
    for i = (g+1):N
        gen = (i - 1) ÷ g
        lo = (gen - 1) * g
        sire[i] = lo + rand(rng, 1:max(1, g ÷ 20))
        dam[i] = lo + g ÷ 2 + rand(rng, 1:(g ÷ 2))
    end
    DataFrame(sire = sire, dam = dam)
end

# progress goes to stderr, unbuffered, so a run that stalls shows where
stage(N, s) = N >= 100_000 && (println(stderr, "  [N = $N] ", s); flush(stderr))

function run(N; seed = 2026)
    rng = Xoshiro(seed)
    ped = pedigree(N; rng)
    g = N ÷ 10
    rec = collect((g+1):N)
    stage(N, "pedigree done; ainv")
    t_ainv = @elapsed Ai = ainv(ped)
    X = sparse(1:length(rec), rand(rng, 1:50, length(rec)), 1.0, length(rec), 50)
    Z = incidence(rec, N)
    y = randn(rng, length(rec))
    stage(N, "mme")
    t_mme = @elapsed m = mme(y, X, [RandomEffect(Z, Ai, 0.3)]; σ²e = 0.7)
    stage(N, "solve")
    t_solve = @elapsed x = solve_mme(m)
    stage(N, "pev (selected inverse)")
    GC.gc()
    b_pev = @allocated t_pev = @elapsed P = pev(m, 1)
    stage(N, "reliability")
    t_rel = @elapsed reliability(m, 1; Kdiag = nrm_diag(ped))
    # check: PEV of 20 random animals from solved columns
    ids = rand(rng, 1:N, 20)
    stage(N, "pev_candidates, 20 columns")
    t_col = @elapsed Γ = pev_candidates(m, 1, ids)
    err = maximum(abs, diag(Γ) .- vec(P)[ids])
    # one REML evaluation (log-likelihood, gradient, AI) on the sparse back-end
    eff = [RandomEffect(Z, Ai, 0.3)]
    GC.gc()
    stage(N, "one REML iteration")
    b_reml = @allocated t_reml = @elapsed r = reml(y, X, eff; σ²e = 0.7, maxiter = 1, backend = :sparse)
    @printf("%9d %9d %8.2f %8.2f %8.2f %8.2f %9.1f %8.2f %9.1f %10.2e %8.2f\n",
        N, nnz(m.lhs), t_ainv, t_mme, t_solve, t_pev, b_pev / 2^20, t_reml, b_reml / 2^20, err,
        t_col / 20 * N)
    flush(stdout)
end

println("# threads ", Threads.nthreads(), ", BLAS ", BLAS.get_num_threads(), ", Julia ", VERSION,
    ", ", Sys.cpu_info()[1].model)
println("#       N   nnz(C)  ainv(s)   mme(s) solve(s)   pev(s)  pev(MiB)  reml(s) reml(MiB) pev maxerr  columns(s)*")
println("# * time to get all PEV by solving one column per animal, extrapolated from 20 columns")
run(1000)                       # compile
for a in ARGS
    run(parse(Int, a))
end
