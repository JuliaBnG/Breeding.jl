function rnd_ped(nn, ng)
    ng > 255 && error("Number of generations must be less than 256")
    ped = DataFrame(
        id = zero(Int32),
        sire = zeros(Int32, nn),
        dam = zeros(Int32, nn),
        sex = rand(Bool, nn),
        grt = UInt8(0),
    )
    for i = 1:ng
        sex = ped.sex[ped.grt .== i-1]
        pm = random_mate(sex, nn)
        acc = Int32((i - 1) * nn)
        off = DataFrame(
            id = zero(Int32),
            sire = pm[:, 1] .+ acc,
            dam = pm[:, 2] .+ acc,
            sex = rand(Bool, nn),
            grt = i,
        )
        append!(ped, off)
    end
    ped.id = Int32.(1:size(ped, 1))
    return ped
end
