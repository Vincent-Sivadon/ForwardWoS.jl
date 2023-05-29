using Zygote
using ForwardWoS

function gradπ_walk(g,∂Ω,x,y,p)
    while true
        d = ∂Ω(x,y)
        if d<1e-3
            return g(x,y,p)
        end
        # Else, take a random pos on the circle of
        # radius d and center xᵢ,yᵢ
        x,y = ForwardWoS.rand_on_ball(x,y,d)
    end
end
function ∂Ω(x,y)
    vx = abs(x)-1.0 ; vy = abs(y)-1.0 ;
    dx = max(vx,0)  ; dy = max(vy,0)
    - (sqrt(dx*dx + dy*dy) + min(max(vx,vy),0))
end
g(x,y,p) = x^2 + p
# p = 3.0
# u = walk(g,∂Ω,0.,0.,p)
# Zygote.gradient(x->walk(g,∂Ω,0.0,0.0,x), p)

function uwithδπs(u,uo,oids,p,xs,ys,ngrid)
    δπ = Dict{Tuple{Int64,Int64},Float64}()

    # Find u and δπ from this π
    for i ∈ 1:ngrid, j ∈ 1:ngrid
        u[i,j],tmpδπ = withgradient(x->gradπ_walk(g,∂Ω,xs[i],ys[j],x), p)
        if (i,j) ∈ oids
            get!(δπ,(i,j),tmpδπ[1])
        end
    end

    # Compute loss function for each known pixels of uo
    for id ∈ oids
        δπ[id] *= 2*(u[id...] - uo[id...])
    end

    u,δπ
end
ngrid = 5
xs = range(-1,1;length=ngrid)
ys = range(-1,1;length=ngrid)
u  = rand(ngrid,ngrid)
p = 2.0
uo = zeros(ngrid,ngrid)
oids = [(3,1),(3,5)]
uo[3,1] = 10.
uo[3,5] = 10.
u,δπ = uwithδπs(u,uo,oids,p,xs,ys,ngrid)

δπ

using CairoMakie
ForwardWoS.plot2D([xs,ys],u)