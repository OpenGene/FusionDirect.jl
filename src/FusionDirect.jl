module FusionDirect

using OpenGene
using OpenGene.Algorithm
using OpenGene.Reference

# package code goes here

# make it compatible for different version of Julia
include("compat.jl")

export detect

import Base: -,
    abs

include("index/index.jl")
include("detect/detect.jl")

end # module
