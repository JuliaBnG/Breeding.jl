module Breeding

using BnGStructs
using DataFrames
using Distributions
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

export sampleID, sampleLoci, eqtl, tbv!, phenotype!, mtpblup
export gene_drop, sum_map, sumMap

end # module Breeding
