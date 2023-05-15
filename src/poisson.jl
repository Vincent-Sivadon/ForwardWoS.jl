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

"Return the one-walk estimator for the poisson problem `p` for position (`x`,`y`)"
function walk(p::PoissonProblem,coords,odim,ϵ)
    uₛ = (odim==1) ? 0. : zeros(odim)  # Source contribution
    while true
        # Closest point to borders
        d = p.∂Ω(coords...)

        # If d<0 then out of bounds
        d < 0 && return (odim==1) ? 0. : zeros(odim)

        # If Distance small enough : return value
        d<ϵ && return p.g(coords...) + uₛ

        # IMPORTANCE SAMPLING
        # -------------------
        # Random point in circle
        incoords = rand_in_ball_is(coords...,d)
        # Add source term
        uₛ += p.f(incoords...) * greens(coords...,incoords...,d) / rpG(sqrt(incoords[1]^2+incoords[2]^2),d)

        # UNIFORM SAMPLING
        # -------------------
        # incoords = rand_in_ball(coords...,d)
        # uₛ += p.f(incoords...) * greens(coords...,incoords...,d) * π*d^2

        # Else, take a random pos on the circle of
        # radius d and center xᵢ,yᵢ
        coords = rand_on_ball(coords...,d)
    end
end