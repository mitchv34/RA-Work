#==========================================================================================
Title: 
Author: Mitchell Valdes-Bobes @mitchv34
Date: 2022-10-07


Description:
==========================================================================================#
# Load packages
using CairoMakie
using LaTeXStrings
include("./model.jl")
#==========================================================================================
plotMatchingSet: Plot the matching set
==========================================================================================#
function plotMatchingSet(prim::Primitives, res::Results; save_path::String="")
    # Unpack necessary objects
    @unpack x_space, y_space, x_lower, x_upper, y_lower, y_upper = prim
    @unpack M = res
    # Create the figure
    fig = Figure(resolution = (500, 500),  font = "CMU Serif", backgroundcolor = :transparent)
    # Create the axis
    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"y", backgroundcolor = :transparent,
                leftspinevisible = true,
                rightspinevisible = false,
                bottomspinevisible = true,
                topspinevisible = false,
                spinewidth = 2
                )
    # Plot the matching set
    hmap = CairoMakie.heatmap!(x_space, y_space, (1 .- M); colormap = :bamako)
    lines!([x_lower, x_upper], [y_upper, y_upper] , color = :black, linewidth = 2, linestyle = :dash)
    lines!([x_upper, x_upper], [y_lower, y_upper] , color = :black, linewidth = 2, linestyle = :dash)
    CairoMakie.xlims!(ax, x_lower, x_upper*1.05)
    CairoMakie.ylims!(ax, y_lower, y_upper*1.05)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    colgap!(fig.layout, 7)

    # Save the figure
    if save_path != ""
        save(save_path, fig)
    end

    return fig # Return the figure
end # function plotMatchingSet
