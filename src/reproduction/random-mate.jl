"""
    random_mate(n1::T, n2::T) where T <: Integer
Randomly selects parents from a population of size `n1` to produce `n2` offspring.
The function returns a tuple containing the parent matrix `pm` and a vector `sex`
indicating the sex of each parent.
"""
function random_mate(n1::T, n2::T) where {T<:Integer}
    all((n1, n2) .> 0) || error("n1 and n2 must be positive integers")
    sex = rand(Bool, n1)
    sirs = rand((1:n1)[sex], n2)
    dams = rand((1:n1)[.!sex], n2)
    pm = sortslices([sirs dams], dims = 1)
    return pm, sex
end

"""
    random_mate(sex::Vector{Bool}, n2::T) where T <: Integer
Randomly selects parents from a population defined by the `sex` vector to produce
`n2` offspring. The function returns a parent matrix `pm` where each row
contains the indices of a sire and a dam.
"""
function random_mate(sex::Vector{Bool}, n2::T) where {T<:Integer}
    n2 > 0 || error("n2 must be a positive integer")
    sirs = rand(findall(sex), n2)
    dams = rand(findall(.!sex), n2)
    pm = sortslices([sirs dams], dims = 1)
    return pm
end
