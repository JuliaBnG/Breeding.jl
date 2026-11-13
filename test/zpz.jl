@testset "Z'Z, batches and haplotypes" begin
    rng = Xoshiro(7)
    n, m = 150, 200 # several 64-bit words and 128-locus tiles
    H = bitrand(rng, m, 2n) # loci × haplotypes
    D = Matrix{Int8}((H[:, 1:2:end] .+ H[:, 2:2:end])') # n × m dosages
    ref = Float64.(D)' * Float64.(D)
    @test zpz_hap!(zeros(m, m), H) == ref
    @test zpz_hap!(zeros(m, m), Matrix(H)) == ref
    @test zpz_hap!(zeros(m, m), Int8.(H)) == ref
    @test_throws ArgumentError zpz_hap!(zeros(m, m), Int8.(2H))
    Hc = falses(n, 2m) # individuals × locus haplotypes
    for k = 1:n, i = 1:m
        Hc[k, 2i-1], Hc[k, 2i] = H[i, 2k-1], H[i, 2k]
    end
    @test zpz_hap!(zeros(m, m), Hc; loci = :cols) == ref
    # BnGStructs containers; m = 200 and n = 150 exercise the 64-bit padding
    @test zpz_hap!(zeros(m, m), Haplotype(H)) == ref
    @test zpz_hap!(zeros(m, m), Genotype(Hc)) == ref
    ZZ = zpz_hap!(zeros(m, m), H[:, 1:100])
    @test zpz_hap!(ZZ, H[:, 101:end]; accumulate = true) == ref
    @test zpz!(zeros(m, m), D; rows = 64) == ref
    @test zpz!(zeros(m, m), Float64.(D)) == ref
    ZZ = zpz!(zeros(m, m), D[1:70, :])
    @test zpz!(ZZ, D[71:end, :]; accumulate = true) == ref
    @test zpz!(zeros(Int, m, m), D) == ref
    @test_throws DimensionMismatch zpz!(zeros(m - 1, m - 1), D)
    @test_throws ArgumentError zpz!(zeros(m, m), D; rows = 0)
    @test_throws ArgumentError zpz_hap!(zeros(m, m), H[:, 1:3])
    @test_throws ArgumentError zpz_hap!(zeros(m, m), H; block = 0)
end
