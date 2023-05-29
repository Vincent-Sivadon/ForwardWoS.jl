module ForwardWoS

using Revise
using Requires
using LinearAlgebra
using CUDA
using CairoMakie
using PrecompileTools

abstract type Problem end

include("utils.jl")
include("solve.jl")
include("laplace.jl")
include("poisson.jl")
include("gpu.jl")
include("withgrad.jl")
include("plotutils/plot2D.jl")

# Function that dispatch on the right type of problem
Problem(g,Ω,∂Ω)   = LaplaceProblem(g,Ω,∂Ω)
Problem(f,g,Ω,∂Ω) = PoissonProblem(f,g,Ω,∂Ω)

# Precompile
# @setup_workload begin
#     # Problem definition
#     g(x,y) = 1. + 2*sin(atan(y,x))
#     ∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
#     Ω = ((-1,1),(-1,1))
#     p = Problem(g,Ω,∂Ω)
#     @compile_workload begin
#         axs,u = SolveGPU(p,ngrid=10,nwalks=50)
#         axs,u = Solve(p,ngrid=10,nwalks=50)        
#         plot2D(axs,u)
#     end

#     # Problem Defintion
#     f(x,y) = 1.  # Uniform heating
#     function sdf(x,y)
#         vx = abs(x)-1 ; vy = abs(y)-1 ;
#         dx = max(vx,0)  ; dy = max(vy,0)
#         - (sqrt(dx*dx + dy*dy) + min(max(vx,vy),0))
#     end
#     Ω = ((-1,1),(-1,1))
#     p = Problem(f,g,Ω,sdf)
#     @compile_workload begin
#         axs,u = SolveGPU(p,ngrid=10,nwalks=50)
#         axs,u = Solve(p,ngrid=10,nwalks=50)
#         plot2D(axs,u)
#     end
# end


# Code loading
# function __init__()
#     @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0"  include("plotutils/plot2D.jl")
#     @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a"     include("plotutils/plot3D.jl")
# end

end # module ForwardWoS
