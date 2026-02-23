"""
    dummy_traits(n::Int, nlc::Int)
Create dummy traits for testing purposes. The function generates `n` dummy
traits with heritability in the range (0.05, 0.55), and σₐ ∈ (0.5, 1.5). A
litter size threshold trait with h² = 0.2 was also added. A random correlation
matrix is generated for these traits. The function returns a tuple containing
the list of traits, the correlation matrix, and the weight of traits that are
included in the total merit index. `nlc` here is to calculate σₐ² from trait
genetic variance. σₐ² = 2σᵤ²/nlc.
"""
function dummy_traits(n::Int, nlc::Int)
    @info "Creating $n dummy traits of h² in (0.05, 0.55)"
    n < 1 || n > 5 && error("n must be in [1, 5]")

    # Dummy traits
    trts, name = [], 'a'  # Dummy trait name
    for _ ∈ 1:n
        t = Trait("t$name"; h² = rand() * 0.5 + 0.05, sex = rand(0:2), σₐ = rand() + 0.5)
        push!(trts, t)
        name += 1
    end
    pp = [0.05, 0.2, 0.35, 0.4]  # threshold category proportions
    push!(trts, Trait("lts", pp; h² = 0.2, sex = 0))  # Litter size trait

    # correlation matrix
    m = n + 1
    R = rand(m, m) * 2 .- 1 # Random correlation matrix
    R = R * R'  # Make it positive definite
    for i = 1:m
        R[i, i] = 4.0
    end
    R /= 4

    pp = rand(n)
    pp /= sum(pp)
    itts = Dict()
    for i = 1:n
        t = trts[i]
        σᵤ² = t.σₐ^2
        σₐ² = 2σᵤ²/nlc
        σₑ² = σᵤ² / t.h² * (1 - t.h²)
        λ = σₑ² / σₐ²
        itts[t.name] = (w = pp[i], λ = λ, h² = t.h²)
    end

    return trts, R, itts
end
