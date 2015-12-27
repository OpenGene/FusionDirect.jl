module FusionDirect

using OpenGene

# package code goes here

import Base: -

include("index/reference.jl")
include("index/index.jl")
include("detect/fusion.jl")
include("detect/detect.jl")

end # module
