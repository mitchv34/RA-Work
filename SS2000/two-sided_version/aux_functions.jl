#==========================================================================================
Title: Axuiliary Functions (Shimer Smith Model - Two Sided Matching)
Author: Mitchell Valdes-Bobes @mitchv34
Date: 2022-10-06

Description: This file contains the auxiliary functions used in the model.
==========================================================================================#
# Load packages
using YAML
#==========================================================================================
setModelPrimitives: Read the model primitives from a YAML file, return a Primitives struct
==========================================================================================#
function setModelPrimitives(model_path)

    # Read the model primitives from a YAML file
    params = YAML.load_file(model_path)
    # Primitives
    r = params["Primitives"]["r"]
    δ = r
    ρ_scale = params["Primitives"]["ρ_scale"] 
    σ_scale = params["Primitives"]["σ_scale"]
    ρ = ρ_scale * r 
    σ = σ_scale * r
    μ = params["Primitives"]["μ"]
    θ = (ρ * μ) / (r + δ)
    ϕ = (σ * (1 - μ)) / (r + δ)
    b =  @eval (x) -> $(Meta.parse(params["Primitives"]["b"]))
    c =  @eval (y) -> $(Meta.parse(params["Primitives"]["c"]))
    # Production function
    prod_funct_params = params["ProdFunct"]["params"]
    if ~isnothing(prod_funct_params) # If the production function has parameters we need to unpack those
        for (k,v) ∈  prod_funct_params
            eval(Meta.parse("$k = $v")) 
        end
    end
    f =  @eval (x, y) -> $(Meta.parse(params["ProdFunct"]["f"]))
    
    # Type_space
    x_lower = params["Type_space"]["x_lower"]
    x_upper = params["Type_space"]["x_upper"]
    y_lower = params["Type_space"]["y_lower"]
    y_upper = params["Type_space"]["y_upper"]
    n_x = params["Type_space"]["n_x"]
    n_y = params["Type_space"]["n_y"]
    x_space = range(x_lower, stop=x_upper,length=n_x)
    y_space = range(y_lower, stop=y_upper,length=n_y)
    # Type_distribution
    x_dist_name = params["Type_distribution"]["x_dist"]
    y_dist_name = params["Type_distribution"]["y_dist"]
    # Depending on the distribution, set the distribution and the cumulative distribution
    # TODO: Add more distributions
    # TODO: @assert that the file hold all the necessary parameters for the distribution
    if x_dist_name == "Uniform"
        x_dist = Uniform(x_lower, x_upper)
    else 
        if isnothing( params["Type_distribution"]["dist_params"]["x"] )
            x_dist = eval(Meta.parse(x_dist_name*"()"))
        end
    end
    if y_dist_name == "Uniform"
        y_dist = Uniform(y_lower, y_upper)
    else
        if isnothing( params["Type_distribution"]["dist_params"]["y"] )
            y_dist = eval(Meta.parse(y_dist_name*"()"))
        end
    end

    L_x(x) = cdf(x_dist, x) 
    l_x(x) = pdf(x_dist, x)
    L_y(x) = cdf(y_dist, x)
    l_y(x) = pdf(y_dist, x)

    # Create the Primitives struct and return it
    return Primitives( r, δ, ρ, σ, μ, θ, ϕ, f, 
                        x_lower,x_upper, x_lower, x_upper, x_dist, y_dist,
                        L_x,l_x,L_y,l_y,n_x,n_y,x_space,y_space,b,c);
end 
