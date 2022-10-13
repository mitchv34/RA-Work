using Distributions
using Plots
using UnicodePlots
using LinearAlgebra

include("funcitons.jl")

# primitives
r = 1.0
δ = r
ρ = 100r 
θ = ρ / (2*(r + δ))
# f(x,y) = (x + y - 1)^2
# f(x,y) = (x + y)^2
f(x,y) = x * y

# Type space
x_lower = 0.0
x_upper = 1.0

# Type distribution
x_dist = Uniform(x_lower, x_upper)
# x_dist = Normal(0.9, 0.01)
L(x) = cdf(x_dist, x)
l(x) = pdf(x_dist, x)

# Divide the space into n_divisions equally spaced points
n_divisions = 500
type_space = range(x_lower, stop=x_upper,length=n_divisions)

# Pre-calculate the production of each pair of types
F = zeros(n_divisions, n_divisions)
for i ∈ 1:n_divisions, j ∈ 1:n_divisions
        F[i,j] = f(type_space[i], type_space[j])
end

##### ----- #####

## Initial matching set
ℳ₀ = Int.(diagm(ones(n_divisions))); # Start with possitive assortative matching
u₀ = ones(n_divisions); # Start with no unmatched agents
w₀ = zeros(n_divisions); # Initialize the value function  u(x) = 0 ∀ x


# Iterate
ℳ₁ = iterate_M(ℳ₀, u₀, w₀);
Plots.heatmap(type_space, type_space, ℳ₁, color=:viridis, legend = false,
        title="Matching Set", xlabel="x", ylabel="y", size = (600, 600))

