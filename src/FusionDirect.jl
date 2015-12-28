module FusionDirect

using OpenGene
using OpenGene.Algorithm

# package code goes here

export detect

import Base: -

include("index/reference.jl")
include("index/index.jl")
include("detect/fusion.jl")
include("detect/detect.jl")

end # module
