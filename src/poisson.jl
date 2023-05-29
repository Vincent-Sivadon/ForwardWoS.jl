"""
Contains field to define a poisson problem
# Arguments
- `f::Function`: source term
- `g::Function`: value imposed on the domain's borders
- `∂Ω::Function`: domain definition (signed distance function)
"""
struct PoissonProblem <: Problem
    f::Function
    g::Function
    Ω::Tuple{Tuple{Float64,Float64}, Tuple{Float64,Float64}}
    ∂Ω::Function
end

"""
    Solve a ForwardWoS problem (Laplace, Poisson) with generic dimensions
    Returns a vector of the axis ranges `axs` (axs[1] -> xrange...) and the solution `u`
"""
function Solve(p::PoissonProblem; ngrid=50, nwalks=200, ϵ=0.001,sample=rand_in_ball,pdf=(d)->1/(π*d^2))
    idim,odim,axs,u = preprocess(p,ngrid)

    # Indices ex 2D with ngrid = 2
    # indices[:] = (1,1), (2,1), (1,2), (2,2)
    ranges = ntuple(x->1:ngrid, idim)
    indices = Iterators.product(ranges...)

    # Solve for each element of u
    for idx in indices
        # Coords associated with idx
        coords = getcoords(axs,idx,idim)

        # Solve at coords
        u[idx...] = walks(p,coords,nwalks,odim,ϵ,sample,pdf)
    end
    axs,u
end

"Return the solution of the laplace problem `p` for position (`xᵢ`,`yᵢ`) with `N` walks"
function walks(p::PoissonProblem,coords,nwalks,odim,ϵ,sample,pdf)
    uᵢ = odim==1 ? 0.0 : zeros(odim)
    for _=1:nwalks
        uᵢ += walk(p,coords,odim,ϵ,sample,pdf)
    end
    uᵢ /= nwalks
end

"Return the one-walk estimator for the poisson problem `p` for position (`x`,`y`)"
function walk(p::PoissonProblem,coords,odim,ϵ,sample,pdf)
    uₛ = (odim==1) ? 0. : zeros(odim)  # Source contribution
    while true
        # Closest point to borders
        d = p.∂Ω(coords...)

        # If d<0 then out of bounds
        d < 0 && return (odim==1) ? 0. : zeros(odim)

        # If Distance small enough : return value
        d<ϵ && return p.g(coords...) + uₛ

        incoords = sample(coords...,d)
        if !isnan(incoords[1])
            pf = pdf(d)
            # rincoords = sqrt((incoords[1]-coords[1])^2 + (incoords[2]-coords[2])^2)
            # pf = rpG(rincoords,d)
            # if pf >= 1e-3
                uₛ += p.f(incoords...) * greens(coords...,incoords...,d) / pf
            # end
        end

        # Else, take a random pos on the circle of
        # radius d and center xᵢ,yᵢ
        coords = rand_on_ball(coords...,d)
    end
end