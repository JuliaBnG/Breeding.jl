@testset "Covariance functions" begin
    G = [132.3 127.0 136.6; 127.0 172.8 200.8; 136.6 200.8 288.0]
    f = covariance_function(G, [90, 160, 240])
    @test f.Ĝ ≈ G
    @test f.M * f.T * f.M' ≈ G
    r = covariance_function(G, [90, 160, 240]; k = 2)
    @test size(r.C) == (2, 2) && r.df == 3
end
