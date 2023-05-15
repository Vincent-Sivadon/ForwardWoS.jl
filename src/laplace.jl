"""
Contains field to define a laplace problem
# Arguments
- `g::Function`: value imposed on the domain's borders
- `∂Ω::Function`: domain definition (signed distance function)
"""
struct LaplaceProblem <: Problem
    g::Function
    Ω::Tuple
    ∂Ω::Function
end

"Return the one-walk estimator for the laplace problem `p` for position (`xᵢ`,`yᵢ`)"
function walk(p::LaplaceProblem,coords,odim,ϵ)
    while true
        # Closest point to borders
        d = p.∂Ω(coords...)

        # If d < 0 then out of bounds
        d < 0 && return (odim==1) ? 0. : zeros(odim)

        # If Distance small enough : return value
        d<ϵ && return p.g(coords...)

        # Else, take a random pos on the circle of
        # radius d and center xᵢ,yᵢ
        coords = rand_on_ball(coords...,d)
    end
end

"Return the one-walk estimator for the laplace problem `p` for position (`xᵢ`,`yᵢ`) using a gradient control-variate strategy"
function walk_withgrad(p::LaplaceProblem,coords,odim,ϵ)
    i = 0 ; coords₁ = coords ; R₁ = 0.
    while true
        # Closest point to borders
        d = p.∂Ω(coords...)

        # If d < 0 then out of bounds
        d < 0 && return (odim==1) ? (0.,coords₁,R₁) : (zeros(odim),coords₁,R₁)

        # If Distance small enough : return value
        d<ϵ && return p.g(coords...),coords₁,R₁

        # Else, take a random pos on the circle of
        # radius d and center xᵢ,yᵢ
        coords = rand_on_ball(coords...,d)

        i==0 && (coords₁=coords ; R₁=d ; i+=1 )
    end
end