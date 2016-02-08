module FusionDirect

using OpenGene
using OpenGene.Algorithm

# package code goes here

export detect

import Base: -,
    abs

include("index/index.jl")
include("detect/detect.jl")

end # module
