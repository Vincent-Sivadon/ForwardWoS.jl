using ForwardWoS
using LinearAlgebra
using Distributions
using Zygote

# Source at (x,y) for a wire of radius rw carrying a current I 
function f(x,y,R=7.0,rw=0.1,I=1.0)
    # sdf (distance to wire)
    d = abs(R - sqrt(x^2+y^2)) 
    if d<=rw
        c = -I/(π*rw^2)
        θ = atan(y,x)
        return c .* [-sin(θ),cos(θ)] # -μ0*j
    end
    return [0.,0.]
end
function sample_f(x,y,d,R=7.,rw=0.1)
    d_to_wire = abs(R - sqrt(x^2 + y^2))
    if d_to_wire < d
        θ = atan(y,x)
        ρ = R + rw*rand()
        rand() < 0.5 && (ρ *= -1)
        return ρ*cos(θ), ρ*sin(θ)
    end
    return NaN,NaN
end
g(x,y) = [0.;0.]
sdf(x,y) = 20.0 - sqrt(x^2 + y^2)
R = 7.

Φᵦ,πs = withgradient(x->abs(-1200-Flux(f,g,sdf,x,100,100)), R)

function Flux(f,g,sdf,R,nsamples=200,nwalks=200)
    Φᵦ = 0.
    for _ ∈ 1:nsamples
        x,y = ForwardWoS.rand_on_ball(0.,0.,R)
        θ = atan(y,x)
        dl = [-sin(θ);cos(θ)]
        A = MagneticPotentialVector(x,y,f,g,sdf,R,nwalks)
        Φᵦ += dot(A,dl)
    end
    Φᵦ /= (nsamples / (2*π*R))
end
function MagneticPotentialVector(x,y,f,g,sdf,R,nwalks)
    A = [0.;0.]
    for _ ∈ 1:nwalks
        A += walk(x,y,f,g,sdf,R)
    end
    A /= nwalks
    return A
end
function walk(x,y,f,g,sdf,R)
    Aₛ = [0.;0.]
    while true
        d = sdf(x,y)
        d<1e-2 && return g(x,y) + Aₛ

        incoords = sample_f(x,y,d,R)
        if !isnan(incoords[1])
            Aₛ += f(incoords...,R) * ForwardWoS.greens(x,y,incoords...,d)
        end
        
        x,y = ForwardWoS.rand_on_ball(x,y,d)
    end
end



pdf(d) = 1.
Ω = ((-20,20),(-20,20))
p = ForwardWoS.Problem(f,g,Ω,sdf)
axs,u = ForwardWoS.Solve(p,ngrid=50,nwalks=10,sample=sample_f,pdf=pdf)
ForwardWoS.plot2D(axs,u,figure=(;resolution=(800,800)))