using ForwardWoS
using CairoMakie

# Problem definition
g(x,y) = 1. + 2*sin(atan(y,x))
∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
Ω = ((-1,1),(-1,1))
p = ForwardWoS.Problem(g,Ω,∂Ω)

# Solving
axs,u = ForwardWoS.SolveGPU(p,ngrid=50,nwalks=200)

# Display
ForwardWoS.plot2D(axs,u)