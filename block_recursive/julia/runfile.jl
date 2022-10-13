
include("model.jl")

prim, res, aux = initializeModel("code/block_recursive/julia/parameters/parameters.yaml")


iterateValueMatch(prim, res)
