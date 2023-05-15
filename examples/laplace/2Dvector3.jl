using ForwardWoS

# Problem definition
g(x,y) = [x;y;x+y]
∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
Ω = ((-1,1),(-1,1))
p = ForwardWoS.Problem(g,Ω,∂Ω)

# Solving
axs,u,∇u = ForwardWoS.Solve_withgrad(p,ngrid=10,nwalks=200)