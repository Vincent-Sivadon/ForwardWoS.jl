module ForwardWoS

using Revise
using Requires
using LinearAlgebra
using CUDA

abstract type Problem end

include("utils.jl")
include("solve.jl")
include("laplace.jl")
include("poisson.jl")
include("gpu.jl")
include("withgrad.jl")

# Function that dispatch on the right type of problem
Problem(g,Ω,∂Ω)   = LaplaceProblem(g,Ω,∂Ω)
Problem(f,g,Ω,∂Ω) = PoissonProblem(f,g,Ω,∂Ω)

# Code loading
function __init__()
    @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0"  include("plotutils/plot2D.jl")
    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a"     include("plotutils/plot3D.jl")
end

end # module ForwardWoS
