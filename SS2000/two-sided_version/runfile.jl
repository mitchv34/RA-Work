#==========================================================================================
Title: 
Author: Mitchell Valdes-Bobes @mitchv34
Date: 2022-10-06


Description:
==========================================================================================#

include("./model.jl")
include("simulation.jl")
include("./ploting_functions.jl")

# for i ∈ 1:4
#     model_path = "./code/SS2000/two-sided_version/params/params_$(i).yaml"
#     println(model_path)
#     prim, res= intializeModel(model_path);
#     iterateMatchingSet(prim, res)
#     plotMatchingSet(prim, res, save_path = "./Notes/(2000) Shimer, Smith/figures/matching_set_$(i).pdf")
# end

i=7
model_path = "./code/SS2000/two-sided_version/params/params_$(i).yaml"
prim, res= intializeModel(model_path);
iterateMatchingSet(prim, res)
plotMatchingSet(prim, res)
 Plots.plot( prim.x_space, res.u )
Plots.plot!( prim.y_space, res.v )
 Plots.plot( prim.x_space, res.U )
Plots.plot!( prim.y_space, res.V )


calculateWages(prim, res)

N_agents  = 60000
agents = obtainInitialAgents(N_agents, prim, res)
