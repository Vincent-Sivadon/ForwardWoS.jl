using Zygote
using CairoMakie
using ForwardWoS

findnearest(ss,v) = findmin(abs.(ss.-v))[2]

function walks(f,g,∂Ω,i,j,nwalks,xs,ys,ngrid,πs)
    u = 0
    for _ ∈ 1:nwalks
        ul = 0. ; us = 0.
        x = xs[i] ; y = ys[j]
        while true
            d = ∂Ω(x,y)
            if d<1e-3
                ul += g(x,y) + us
                break
            end

            inx,iny = ForwardWoS.rand_in_ball(x,y,d)
            # Closest index from x,y
            i_fromxy = findnearest(ys,iny)
            j_fromxy = findnearest(xs,inx)
            us += f(inx,iny,πs[(i_fromxy-1)*ngrid + j_fromxy]) * ForwardWoS.greens(x,y,inx,iny,d) * π*d^2

            # Else, take a random pos on the circle of
            # radius d and center xᵢ,yᵢ
            x,y = ForwardWoS.rand_on_ball(x,y,d)
        end
        u += ul
    end
    u /= nwalks
    return u
end

function InverseSolve(uo,f,ngrid,nwalks)
    ∇π = zeros(ngrid*ngrid)
    πs = rand(ngrid*ngrid)
    u = zeros(ngrid,ngrid)
    xs = range(-1,1;length=ngrid)
    ys = range(-1,1;length=ngrid)

    # while max(abs.(∇π)...) > 1e-1
    for iter ∈ 1:500
        ∇π .= 0.
        for (idx,uoi) ∈ uo
            u[idx...],tmp∇π = withgradient(x->walks(f,g,∂Ω,idx[1],idx[2],nwalks,xs,ys,ngrid,x), πs)
            ∇π .+= ((2*(uoi - u[idx...])) .* tmp∇π[1])
        end
        πs .+= 0.5 .* ∇π

        if iter % 50 == 0
            axs,us = Solve(πs,f,g,∂Ω,ngrid,nwalks)
            # display(ForwardWoS.plot2D([xs,ys],reshape(πs,(ngrid,ngrid))))
            display(ForwardWoS.plot2D([xs,ys],us))
        end
    end
    πs
end

function Solve(πs,f,g,∂Ω,ngrid,nwalks)
    u = zeros(ngrid,ngrid)
    xs = range(-1,1;length=ngrid)
    ys = range(-1,1;length=ngrid)
    for i ∈ 1:ngrid, j ∈ 1:ngrid
        u[i,j] = walks(f,g,∂Ω,i,j,nwalks,xs,ys,ngrid,πs)
    end
    [xs,ys],u
end

function ∂Ω(x,y)
    vx = abs(x)-1.0 ; vy = abs(y)-1.0 ;
    dx = max(vx,0)  ; dy = max(vy,0)
    - (sqrt(dx*dx + dy*dy) + min(max(vx,vy),0))
end
g(x,y) = 0.
fs(x,y,p) = p
nwalks = 100
ngrid = 20
uo = Dict(
    (10,10)=>10.,
    (10,11)=>10.,
    (11,10)=>10.,
    (11,11)=>10.
)

πs = InverseSolve(uo,fs,ngrid,nwalks)

πs = ones(ngrid,ngrid)
axs,us = Solve(πs,fs,g,∂Ω,ngrid,nwalks)
fig = ForwardWoS.plot2D(axs,us)