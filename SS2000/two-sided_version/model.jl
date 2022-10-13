#==========================================================================================
Title: Shimer Smith Model - Two Sided Matching
Author: Mitchell Valdes-Bobes @mitchv34
Date: 2022-10-05

Description: This file contains the necessary functions to compute the equilibrium of the
                two-sided labor market verison of the Shimer Smith model.
==========================================================================================#
# Load packages
using Distributions
using Plots
using LinearAlgebra
using Parameters
include("./aux_functions.jl") # Auxiliary functions
#==========================================================================================
Model: Structure holding the preimitives od the model. 
==========================================================================================#
@with_kw struct Primitives
    # Primitives
    # TODO: Add parameter descriptions
    r::Float64  #
    δ::Float64  #
    ρ::Float64  #
    σ::Float64  #
    μ::Float64  #
    θ::Float64  #
    ϕ::Float64  #
    f::Function #

    # Type space
    x_lower::Float64 #
    x_upper::Float64 #
    y_lower::Float64 #
    y_upper::Float64 #

    # # Type distribution
    x_dist#::Distribution = Uniform(x_lower, x_upper)
    y_dist#::Distribution = Uniform(y_lower, y_upper)
    L_x::Function #
    l_x::Function #
    L_y::Function #
    l_y::Function #

    # # Divide the space into n_divisions equally spaced points
    n_x::Int64
    n_y::Int64
    x_space::Vector{Float64} = range(x_lower, stop=x_upper,length=n_x)
    y_space::Vector{Float64} = range(y_lower, stop=y_upper,length=n_y)

    b::Function #
    c::Function #
end # struct Primitives
#==========================================================================================
Resutls: Structure holding the results of the model. 
==========================================================================================#
mutable struct Results
    M::Matrix{Float64}                  # Matching Set
    u::Vector{Float64}                  # Unemployed worker distribution
    v::Vector{Float64}                  # Unfilled vacancy distribution
    U::Vector{Float64}                  # Unemployed worker value function
    V::Vector{Float64}                  # Unfilled vacancy value function 
    W::Matrix{Float64}                  # Wages W(x,y) is the wage for a worker with type x and a vacancy with type y.

    # Auxiliary variables useful to precompute
    F::Matrix{Float64}                  # Production of each pair of types

    # Construtor
    function Results(n_x::Int64, n_y::Int64)
        M = zeros(n_x, n_y) 
        u = zeros(n_x)
        v = zeros(n_y)
        U = zeros(n_x)
        V = zeros(n_y)
        W = zeros(n_x, n_y)
        F = zeros(n_x, n_y)
        new(M, u, v, U, V, W, F)
    end # function Results
end # struct Results
#==========================================================================================
intializeModel: Initialize the model. Recoeves a path to the model primitives and returns 
            primitives and results structures.
==========================================================================================#
function intializeModel(path::String)
    # Load primitives
    prim = setModelPrimitives(path)
    # Initialize results
    res = Results(prim.n_x, prim.n_y)
    return prim, res
end # function initialize
#==========================================================================================
updateF!:  Precomputes the production matrix F. F[i,j] = f(x_i, y_j)
==========================================================================================#
function updateF!(prim::Primitives, res::Results)
    @unpack f, n_x, n_y = prim
    F = zeros(prim.n_x, prim.n_y)
    for i ∈ 1:prim.n_x, j ∈ 1:prim.n_y
            F[i,j] = prim.f(prim.x_space[i], prim.y_space[j])
    end
    res.F = F
end # function updateF!
#==========================================================================================
getNewDistributions: Given a mathching set and intial distribution of unemployed workers, 
                    and empty vacancies, returns the new distribution of unemployed workers
                    and empty vacancies.
==========================================================================================#
function getNewDistributions(prim::Primitives, res::Results)
    
    # Unpack primitives
    @unpack n_x, n_y, x_space, y_space, l_x, l_y, δ, σ, ρ = prim
    
    # Initialize the new distributions
    new_u = zeros(n_x)
    new_v = zeros(n_y)

    # Iterate over the x grid
    for i_x ∈ 1:n_x
        # u(x) =  δ lₓ(x) / (δ + σ ∫ M(x,y) u(y) dy)
        new_u[i_x] = δ .* l_x(x_space[i_x]) ./( δ .+ σ .* (sum(res.M[i_x, :] .* res.v) / n_y)  )
    end

    # Iterate over the y grid
    for i_y ∈ 1:n_y
        # v(y) =  δ lᵧ(y) / (δ + σ ∫ M(x,y) v(x) dx)
        new_v[i_y] =  δ .* l_y(y_space[i_y]) ./( δ .+ ρ .* (sum(res.M[:, i_y] .* res.u) / n_x)  )
    end

    return new_u, new_v

end # function getNewDistributions
#==========================================================================================
iterateDistributions: Given a mathching set and intial distribution of unemployed workers, 
                    iterates until the distributions converge.
==========================================================================================#
function iterateDistributions(prim::Primitives, res::Results;
                                n_iter_max = 10000,  tol::Float64 = 1e-10)
    residual = Inf; # Initialize the residual to Inf
    n_iter = 0; # Initialize the number of iterations to 0
    while (residual > tol) && (n_iter < n_iter_max)
        n_iter += 1; # Update the number of iterations
        u_next, v_next = getNewDistributions(prim, res) # Get the new distributions
        # calculate the residual
        residual_u = norm(u_next[:] .- res.u[:])
        residual_v = norm(v_next[:] .- res.v[:])
        residual = max(residual_u, residual_v)
        # @show residual
        # Update the results
        res.u = copy( u_next )
        res.v = copy( v_next )
    end
end # function iterateDistributions
#==========================================================================================
getNewValueFunctions: Given a mathching set, distribution of unemployed workers, 
                    and empty vacancies, and initial gues of value functions returns the 
                    new value functions. 
==========================================================================================#
function getNewValueFunctions( prim::Primitives, res::Results )
    
    # Unpack primitives
    @unpack n_x, n_y, x_space, y_space, δ, σ, ρ, θ, ϕ, f, b, c = prim
    # Unpack results
    @unpack M, u, v, F = res
    # Initialize the new value functions
    new_V = zeros(n_x)
    new_U = zeros(n_y)
    # Iterate over the x grid
    for i_x ∈ 1:n_x
        den = 1 .+ θ .* sum(res.M[i_x, :] .* v) / n_y
        num = b(x_space[i_x]) .+ θ .* sum( (M[i_x, :] .* (F[i_x, :]  .- res.V)) .* v ) / n_y
        new_U[i_x] = num / den  
    end
    # Iterate over the y grid
    for i_y ∈ 1:n_y
        den = 1 + ϕ .* sum( M[:, i_y] .* u) / n_x
        num = -c(y_space[i_y]) .+ ϕ .* sum( (M[:, i_y] .*(F[:, i_y]  .- res.U)) .* u ) / n_x
        new_V[i_y] = num / den
    end
    
    return new_U, new_V

end # function getNewValueFunctions
#==========================================================================================
iterateValueFunctions: Given a mathching set, distribution of unemployed workers, 
                    and empty vacancies, and initial gues of value functions iterates until 
                    the value functions converge. 
==========================================================================================#
function iterateValueFunctions(prim::Primitives, res::Results; n_iter_max = 10000,  tol::Float64 = 1e-10)
    residual = Inf; # Initialize the residual to Inf
    n_iter = 0; # Initialize the number of iterations to 0
    while (residual > tol) && (n_iter < n_iter_max)
        n_iter += 1; # Update the number of iterations
        U_next, V_next = getNewValueFunctions(prim, res) # Get the new value functions
        # calculate the residußal
        residual_U = norm(U_next[:] .- res.U[:])
        residual_V = norm(V_next[:] .- res.V[:])
        residual = max(residual_U, residual_V)
        # @show residual
        res.U = copy( U_next )
        res.V = copy( V_next )
    end
    # Update the results

end # function iterateValueFunctions
#==========================================================================================
iterateMatchingSet: Given a mathching set, distribution of unemployed workers, 
                    and empty vacancies, and initial gues of value functions iterates until 
                    the matching set converges.
==========================================================================================#
function iterateMatchingSet( prim::Primitives, res::Results ;n_iter_max = 10000,  tol::Float64 = 1e-10)

    # Update the productivity grid
    updateF!(prim, res)

    residual = Inf;
    n_iter = 0;


    while (residual > tol) && (n_iter < n_iter_max)
        n_iter += 1;
        # Iterate to obtain distributions consistent with the matching set
        iterateDistributions(prim, res)
        # Iterate to obtain value functions implied by the matching set and distributions
        iterateValueFunctions(prim, res)
        # Calculate the new matching set implied by the value functions
        new_M = Int.(res.F .- res.U' .- res.V .≥ 0)
        # Calculate the residual
        residual = norm(new_M[:] .- res.M[:])
        @show residual
        # Update the matching set
        res.M = copy( new_M )
    end

end # function getNewMatchingSet
#==========================================================================================
calculateWages: Given a mathching set, distribution of unemployed workers, 
                and empty vacancies, calculates the wages of any pair (x,y).
==========================================================================================#
function calculateWages(prim::Primitives, res::Results)
    # Make sure that F is up to date
    updateF!(prim, res)
    # Unpack primitives
    @unpack n_x, n_y, μ = prim
    # Unpack results
    @unpack F, U, V = res
    # Initialize the wage matrix
    W = zeros(n_x, n_y)
    for i_x ∈ 1:n_x # Iterate over the x grid
        for i_y ∈ 1:n_y # Iterate over the y grid
            # W(x,y) =  μ[f(x,y)-U(x)-V(x)] + U(x) 
            W[i_x, i_y] = μ * (F[i_x, i_y] - U[i_x] - V[i_y]) + U[i_x]
        end
    end
    # Update the results
    res.W = copy( W )
end # function calculateWages