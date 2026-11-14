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
