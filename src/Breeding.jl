module Breeding

using BnGStructs
using DataFrames
using LinearAlgebra
using Statistics
using Random
using SparseArrays
using StatsBase

include("util.jl")
include("reproduction/reproduce.jl")
include("prediction/predict.jl")
include("trait/trait.jl")
include("ocs/ocs.jl")

end # module Breeding
