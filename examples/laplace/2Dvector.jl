using ForwardWoS

# Problem definition
g(x,y) = [x;y]
∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
Ω = ((-1,1),(-1,1))
p = ForwardWoS.Problem(g,Ω,∂Ω)

# Solving
axs,u = ForwardWoS.Solve(p,ngrid=25,nwalks=200)

# Display
using CairoMakie
ForwardWoS.plot2D(axs,u)