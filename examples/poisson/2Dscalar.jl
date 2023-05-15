using ForwardWoS
using CairoMakie

# Problem Defintion
s(x,y) = 1.  # Uniform heating
g(x,y) = 0.  # Cold borders
function sdf(x,y,Lx=1,Ly=1)
    vx = abs(x)-Lx ; vy = abs(y)-Ly ;
    dx = max(vx,0)  ; dy = max(vy,0)
    - (sqrt(dx*dx + dy*dy) + min(max(vx,vy),0))
end
Ω = ((-1,1),(-1,1))
p = ForwardWoS.Problem(s,g,sdf,Ω)

# Solving
axs,u = ForwardWoS.SolveGPU(p,ngrid=200,nwalks=1000)

# Display
fig, = ForwardWoS.plot2D(axs,u)
fig

# surface(x,y,u)
# Analytical Solution
function ua(x,y,n_modes=1)
    sol = (1-x^2)/2
    for k=1:n_modes
        if (k%2 == 0) continue; end
        sol -= 16/π^3 *
            sin(k*π*(1+x)/2)/(k^3*sinh(k*π)) *
            (sinh(k*π*(1+y)/2) + sinh(k*π*(1-y)/2))
    end
    sol
end
heatmap(axs...,ua,aspect_ratio=1)