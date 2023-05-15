using ForwardWoS
using GLMakie

# Problem definition
g(x,y,z) = [x;y;z]
∂Ω(x,y,z) = 1.0 - sqrt(x*x + y*y + z*z)
Ω = ((-1,1),(-1,1),(-1,1))
p = ForwardWoS.Problem(g,Ω,∂Ω)

# Solving
ngrid=15
axs,u = ForwardWoS.Solve(p,ngrid=ngrid,nwalks=200)

# Display
ForwardWoS.plot3D(axs,u)