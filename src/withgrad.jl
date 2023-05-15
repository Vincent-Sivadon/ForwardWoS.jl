function Solve_withgrad(p::Problem; ngrid=50, nwalks=200, ϵ=0.001)
    idim,odim,axs,u,∇u = preprocess_withgrad(p,ngrid)

    # Indices ex 2D with ngrid = 2
    # indices[:] = (1,1), (2,1), (1,2), (2,2)
    ranges = ntuple(x->1:ngrid, idim)
    indices = Iterators.product(ranges...)

    # Solve for each element of u
    for idx in indices
        # Coords associated with idx
        coords = getcoords(axs,idx,idim)

        # Solve at coords
        u[idx...],∇u[idx...] = walks_withgrad(p,∇u,coords,nwalks,odim,idim,ϵ)
    end
    
    axs,u,∇u
end

function walks_withgrad(p::Problem,∇u,coords,nwalks,odim,idim,ϵ)
    uᵢ = odim==1 ? 0.0 : zeros(odim)
    ∇uᵢ = zeros(size(∇u[1]))
    for j=1:nwalks
        # j-th walk informations
        uᵢⱼ,coords₁,R₁ = walk_withgrad(p,coords,odim,ϵ)

        # Add this walk contribution to the solution
        uᵢ += uᵢⱼ ./ nwalks

        # Gradient
        R₁ <= 100*ϵ && continue
        diff = vec([(coords₁.-coords)...])
        curr∇uᵢ = (2 .* (uᵢⱼ.-uᵢ)./R₁^2) * diff' ./ nwalks
        odim==1 && (curr∇uᵢ = curr∇uᵢ' ; curr∇uᵢ=reshape(curr∇uᵢ,(idim,)))
        ∇uᵢ += curr∇uᵢ
    end
    uᵢ,∇uᵢ
end