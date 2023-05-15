# WalkOnSphereSim.jl
This package provides tools to solve PDE (currently Laplace and Poisson problems) with the Monter Carlo Walk On Sphere strategy.

# Example
![OneWalk](Figures/OneWalk.svg)

# Usage
Here is the core of how to use this package to solve PDEs.
We will in this example solve the following Poisson problem :

$$\Delta u(x,y) = f(x,y) \quad \forall (x,y)\in \Omega \qquad \text{and} \qquad u(x,y) = g(x,y) \quad \forall (x,y) \in \partial \Omega$$

```julia
using ForwardWoS
using CairoMakie

# Problem definition
g(x,y) = 1. + 2*sin(atan(y,x))
∂Ω(x,y) = 1.0 - sqrt(x*x + y*y)
Ω = ((-1,1),(-1,1))
p = ForwardWoS.Problem(g,∂Ω,Ω)

# Solving
axs,u = ForwardWoS.SolveGPU(p,ngrid=50,nwalks=200)

# Display
ForwardWoS.plot2D(axs,u)
```

The domain's borders ($\partial \Omega$) is defined by a signed distance function.
Because we gavee the `ForwardWoS.Problem` a function f, the solver will use a strategy considering it's a Poisson problem. Calling `ForwardWoS.Problem` with 3 arguments will imply a Laplacian problem. The rest stays the same.

# Generic Dimensions

The dimensions are handled in a generic way, so the possible modelisations are : 2D grid & scalar solution, 3D grid & 2D vector solution... The GPU version can for now only handle 2D grid & scalar solution. Some examples below (the plots are also generic on the solution dimension, but you must specify ForwardWoS.plot2D or plot3D because it will then use CairoMakie or GLMakie).

![OneWalk](Figures/2Dscalar.png)
![OneWalk](Figures/2Dvector.png)
![OneWalk](Figures/3Dscalar.png)
![OneWalk](Figures/3Dvector.png)

# Bibliography
Mainly based on the paper :
Rohan Sawhney and Keenan Crane. 2020. Monte Carlo Geometry Processing:
A Grid-Free Approach to PDE-Based Methods on Volumetric Domains. ACM
Trans. Graph. 38, 4, Article 1 (July 2020), 18 pages