using ForwardWoS
using CairoMakie

# Problem definition
g(x,y) = 1. + 2*sin(atan(y,x))
∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
Ω = ((-1,1),(-1,1))
p = ForwardWoS.Problem(g,∂Ω,Ω)

# Solving
axs,u,∇u = ForwardWoS.Solve_withgrad(p,ngrid=20,nwalks=10)

# Display
f, = ForwardWoS.plot2D(axs,∇u)
f