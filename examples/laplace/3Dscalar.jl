using ForwardWoS
using GLMakie

# Problem definition
g(x,y,z) = z + 2*sin(atan(y,x))
∂Ω(x,y,z) = 1.0 - sqrt(x*x + y*y + z*z)
Ω = ((-1,1),(-1,1),(-1,1))
p = ForwardWoS.Problem(g,Ω,∂Ω)

# Solving
axs,u = ForwardWoS.Solve(p,ngrid=25,nwalks=200)

# Display
ForwardWoS.plot3D(axs,u)