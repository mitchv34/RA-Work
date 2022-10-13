#==========================================================================================
Title:  Model
Author: Mitchell Valdes-Bobes @mitchv34
Date:   29 September 2022

Description: 
==========================================================================================#

# Import packages
using Parameters
using LinearAlgebra
using YAML
using Term # For fancy printing, you don't need this

# Import auxiliary functions
include("tauchen.jl")


#==========================================================================================
Primitives: Structure tha holds the primitives of the model
==========================================================================================#
@with_kw struct Primitives
    # Parameters
    β      ::Float64 = 0.996    # discount factor
    ρ      ::Float64 = 0.98     # persistence of productivity
    σₑ     ::Float64 = 0.01     # standard deviation of productivity shocks
    n_Z    ::Int64   = 15       # number of productivity states
    γ      ::Float64 = 0.6      # job finding rate parameter
    κ      ::Float64 = 2.37     # vacancy posting cost
    δ_bar  ::Float64 = 0.012    # job destruction rate
    ## utility function parameters:
    σ      ::Float64 = 0.0 
    α      ::Float64 = 1.0
    χ      ::Float64 = 2.0
    a      ::Float64 = 1/3       # job posting opportunities parameter
    b_prob ::Float64 = 0.1       # probability of loosing unemployment benefits
    # Grids
    n_S    ::Int64               # size of search intensity grid
    s_min  ::Float64             # minimum search intensity
    s_max  ::Float64             # maximum search intensity
    s_grid ::Vector{Float64}    # search intensity grid
    n_W    ::Int64               # number of wage grid points
    W_max  ::Float64             # maximum wage
    W_min  ::Float64             # minimum wage
    W_grid ::Vector{Float64}     # wage grid
    Z_grid ::Vector{Float64}     # productivity grid
    δ_grid ::Matrix{Float64}     # job destruction grid
    Π      ::Matrix{Float64}     # transition matrix for productivity shocks
    

    # Functions 
    p      ::Function
    λ      ::Function
    u      ::Function
    
end # struct Primitives
    
#==========================================================================================
getPrimitives: Reads primitives from a file and returns a Primitives object
==========================================================================================#
function getPrimitives(parameters_path::String)


    # Read parameters from YAML file
    parameters = YAML.load_file(parameters_path)["Primitives"]

        
    # Retrieve parameters
    β         = parameters["beta"]
    ρ         = parameters["rho"]
    σₑ        = parameters["sigma_e"]
    n_Z       = parameters["n_z"]
    γ         = parameters["gamma"]
    κ         = parameters["kappa"]
    δ_bar     = parameters["delta_bar"]
    σ         = parameters["sigma"]
    α         = parameters["alpha"]
    χ         = parameters["chi"]
    a         = parameters["a"]
    b_prob    = parameters["b_prob"]
    s_min     = parameters["s_min"]
    s_max     = parameters["s_max"]
    n_s       = parameters["n_s"]
    s_grid    = range(s_min, s_max, length=n_s)
    n_W       = parameters["n_W"]
    Z_grid, Π = tauchenMethod(0.0, σₑ, ρ, n_Z; q = 3)
    Z_grid    = exp.(Z_grid)
    W_max     = maximum(Z_grid) #parameters["W_max"]
    W_min     = minimum(Z_grid)#parameters["W_min"]
    W_grid    = range(W_min, W_max, length=n_W)
    δ_grid    = zeros(n_Z, n_W)

    for i ∈ eachindex(Z_grid), j ∈ eachindex(W_grid)
        δ_grid[i, j] = (Z_grid[i] ≥ W_grid[j]) ?  δ_bar : 1.0
    end

    p  = (θ) -> θ(1 + θ^γ)^(-1/γ) # job finding rate function
    λ  = (s) -> s^a # job posting opportunities function
    u  = (c,s) -> (c^(1-σ) - 1)/(1-σ) - α * s^χ  # utility function``


    # Create Primitives object
    Primitives(
    β, ρ, σₑ, n_Z, γ,
    κ, δ_bar, σ, α, χ, 
    a, b_prob,
    n_s,
    s_min,
    s_max,
    n_s,
    s_grid,
    n_W   ,
    W_max ,
    W_min,
    W_grid,
    Z_grid,
    δ_grid, Π, p, λ, u)

end # function getPrimitives

#==========================================================================================
Results: Structure tha holds the resuts of the model
==========================================================================================#
mutable struct Results
    # Value functions 
    J::Array{Float64, 2} # value of an ongoing match J(w,z)
    θ::Array{Float64, 2} # value of a new match θ(w,z)
    W::Array{Float64, 3} # value of an employedd household W(w,z,μ)
    U::Array{Float64, 3} # value of an unemployed household U(w,z,μ)
    # Policy functions
    s::Array{Float64, 3} # job search intensity s(w,z,μ)
end # struct Results

#==========================================================================================
Aux: Structure tha holds the auxiliary objects of the model that are precomputed
==========================================================================================#
mutable struct Aux
    # Value functions 
    ZW::Array{Float64, 2} # Difference between productivity and wage
    U ::Array{Float64, 2} # Utility of every possible combination of wage and search intensity
end # struct Results



#==========================================================================================
initializeModel: Initializes the model
==========================================================================================#
function initializeModel(parameters_path::String)
    
        # Get primitives
        prim = getPrimitives(parameters_path)
    
        # Initialize results
        ## Value functions
        J = zeros(prim.n_Z, prim.n_W)
        θ = zeros(prim.n_Z, prim.n_W)
        W = zeros(prim.n_Z, prim.n_W, 2)
        U = zeros(prim.n_Z, prim.n_W, 2)
        ## Policy functions
        s = zeros(prim.n_Z, prim.n_W, 2)

        # Precompute auxiliary objects
        ZW = hcat([prim.Z_grid .- prim.W_grid[i] for i ∈ eachindex(prim.W_grid)]...)
        U = zeros(prim.n_W, prim.n_S)
        for i ∈ eachindex(prim.W_grid), j ∈ eachindex(prim.s_grid)
            U[i, j] = prim.u(prim.W_grid[i], prim.s_grid[j])
        end
        return prim, Results(J, θ,  W, U, s), Aux(ZW, U)
    
end # function initializeModel

#==========================================================================================
getNewJ: Apply the Bellman operator to J
==========================================================================================#
function getNewJ(prim::Primitives, res::Results, aux::Aux)
    
    # Unpack primitives
    @unpack W_grid, Z_grid, Π, δ_grid = prim
    # Unpack ZW matrix
    @unpack ZW = aux

    J_new = zeros(prim.n_Z, prim.n_W)
    # Iterate over wage grid
    r = 0.01
    for w_i in 1:n_W
        
        J_new[:, w_i] = ZW[:, w_i]' .+ ((1 .- δ_grid[:, w_i]) .* res.J[:, w_i] )' * Π/(1 + r)
        
    end # for w in 1:n_W
    
    return J_new

end # function getNewJ

#==========================================================================================
iterateValueMatch: Solves the the value fun of an ongoing match
==========================================================================================#
function iterateValueMatch(prim::Primitives, res::Results; n_iter_max::Int64=5000, tol::Float64=1e-6)
    
    println(@bold @blue "Iterating value of an ongoing match...")

    # Unpack primitivesβ, W_grid, Z_grid, Π, δ_grid = prim
    
    # Iterate
    dist = Inf
    n_iter = 0
    while dist > tol && n_iter < n_iter_max
        n_iter += 1
        # Update old value function
        J_new = getNewJ(prim, res)
        
        # Update distance
        dist = norm(J_new - res.J)
        # Update iteration counter
        n_iter += 1
        
        # Update value function
        res.J = copy(J_new)

    end # while dist > tol && n_iter < n_iter_max
    
    # Print results
    if dist < tol
        println(@bold @green "Value function of an ongoing match converged after $n_iter iterations tolerance: $tol")
    else
        println(@bold @red "Value function of an ongoing match did not converge after $n_iter_max iterations tolerance: $tol")
    end # if dist < tol
end # function iterateValueMatch

#==================================================================================================
solveSubMarketTightness: Compute each submarket tightness
==================================================================================================#
function solveSubMarketTightness(prim::Primitives, res::Results)
    
    @unpack κ = prim

    res.θ = κ / res.J 
    res.θ *= Int.(res.J .≥ 0)
    
end # function solveSubMarketTightness

#==========================================================================================
getNewW_U: Apply the Bellman operator to W and U
==========================================================================================#
function getNewW_U(prim::Primitives, res::Results, aux::Aux)
    
    # Unpack primitives
    @unpack n_W, n_Z = prim
    # Unpack U matrix
    @unpack U = aux
    # Unpack results
    @unpack θ, s = res

    
    # Iterate over wage grid
    for w_i ∈ 1:n_W # Iterate over wage grid
        u_emp = U[w_i, 1] # utility of being employed 
        for b_i ∈ 1:n_w # Iterate over wage grid (now as benefit)
            for z_i ∈ 1:n_Z
                
            end
        end
    end # for w in 1:n_W


    return W_new

end # function getNewW
