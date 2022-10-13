#==========================================================================================
Title: Simulaiton
Author: Mitchell Valdes-Bobes @mitchv34
Date: 2022-10-07

Description: This file contains the functions to simulate a two-sided Shimer-Smith model.
==========================================================================================#
# Load packages
using Distributions
using DataFrames
include("./model.jl")
#==========================================================================================
Agent: Sructure that defines an agent.
==========================================================================================#
mutable struct Agent
    # Parameters
    age ::Int64                # Age of the agent (in model periods, 0 at birth)
    X   ::Float64              # Agent's productivity type.
    Y   ::Vector{Float64}      # Employment history ( yₜ = y≥0 if employed at firm y,
                                                    # yₜ = -1 if unemployed)
    W	::Vector{Float64}	   # Wage history (wₜ = w if employed at firm yₜ, wₜ = b(x) if unemployed)

    # Constructor
    function Agent(X::Float64)
        new(0, # Initial age is 0
            X, # Productivity type
            Vector{Float64}(undef, 0), # Empty employment history
            Vector{Float64}(undef, 0) # Empty wage history
        )
    end
end # struct Agent
# TODO: Interpolate the wage function to be able to simulate wages for any (x,y) combination.
# TODO: Interpolat the matching set function to be able to simulate the matching probability for any (x,y) combination.

#==========================================================================================
description: 
==========================================================================================#
function obtainInitialAgents( n_agents::Int, prim::Primitives, res::Results)
    # Unpack necessary objects
    @unpack x_dist = prim
    # Pre allocate space for the agents
    agents = [Agent(rand(x_dist)) for i in 1:n_agents]
    # for i \in 1:n_agents
    #     # Draw the agent's productivity type
    #     agent.X = 
    #     # Draw the agent's initial employment status
    # end 
end # function obtainInitialAgents
#==========================================================================================
getAgentsData: Creates a DataFrame with the agents' data. 
==========================================================================================#
function getAgentsData(agents::Vector{Agent})
    
    # Create vectors to store the data
    
    for i ∈ 1:length(agents)
        # Unpack the agent
        agent = agents[i]
        # Add the agent's data to the DataFrame
        df[i, :] = [agent.age, agent.X, agent.Y, agent.W]
    end
    return df
end # function getAgentsData