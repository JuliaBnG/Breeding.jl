using Breeding
using DataFrames
using RelationshipMatrices
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using Test
using Aqua
using BnGStructs: Trait, Haplotype, Genotype

# A⁻¹ from RelationshipMatrices for parent vectors
ainv_ref(s, d) = Matrix(ainv(DataFrame(sire = s, dam = d)))

# Mrode Example 4.1
const SIRE = [0, 0, 0, 1, 3, 1, 4, 3]
const DAM = [0, 0, 0, 0, 2, 2, 5, 6]
const CALF = [4, 5, 6, 7, 8]
const SEX = ["M", "F", "F", "M", "M"]
const WWG = [4.5, 2.9, 3.9, 3.5, 5.0]

# Shared by the files below: the animal model of Example 4.1
Ai = ainv_ref(SIRE, DAM)
X, lv = incidence(SEX; levels = ["M", "F"])
Z = incidence(CALF, 8)
a = RandomEffect(Z, sparse(Ai), 20.0)

@testset "Breeding" begin
    @testset "Aqua" begin
        Aqua.test_all(Breeding)
    end
    include("ocs.jl")
    include("design.jl")
    include("mme.jl")
    include("multicomponent.jl")
    include("mtpblup.jl")
    include("zpz.jl")
    include("iod.jl")
    include("ram.jl")
    include("pedigree.jl")
    include("genomic.jl")
    include("bayes.jl")
    include("threshold.jl")
    include("covfun.jl")
    include("vc.jl")
    include("reml.jl")
end
