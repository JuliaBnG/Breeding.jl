@testset "Threshold models" begin
    # Mrode Example 15.1: first iteration from the book
    herd = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    sex = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1]
    sire = [1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4]
    N = [
        1 0 0;
        1 0 0;
        1 0 0;
        0 1 0;
        1 0 1;
        3 0 0;
        1 1 0;
        0 1 0;
        1 0 0;
        2 0 0;
        1 0 0;
        0 0 1;
        1 0 1;
        1 0 0;
        0 1 0;
        0 0 1;
        0 1 0;
        1 0 0;
        2 0 0;
        2 0 0
    ]
    Xt = [first(incidence(herd)) first(incidence(sex))]
    Zt = first(incidence(sire))
    r = threshold_model(
        N,
        Xt,
        Zt,
        ainv_ref([0, 0, 1, 3], [0, 0, 0, 0]),
        19;
        t0 = [0.468, 1.080],
        zero = [1, 3],
    )
    @test r.history[:, 1] ≈ [
        0.441008,
        1.044792,
        0,
        0.286869,
        0,
        -0.358323,
        -0.041528,
        0.057853,
        0.039850,
        -0.065178,
    ] atol = 1e-5
    @test round.([r.t; r.u], digits = 4) ≈
          [0.4378, 1.0675, -0.0434, 0.0592, 0.0412, -0.0660]
    P = category_probabilities(r.t, [0.0, 1.0])
    @test sum(P; dims = 2) ≈ ones(2)
end
