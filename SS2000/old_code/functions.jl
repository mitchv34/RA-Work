
function get_new_u(old_u::Vector, M::Matrix)
    # Given a Match set and Starting with initial unmatched rates iterate 
    # Calcualte the inegral in the denominator
    new_u = zeros(n_divisions)
    for i ∈ 1:n_divisions
        den = δ + ρ * sum( M[i,:] .* old_u)
        # Calcualte next unmatched rates
        new_u[i] = δ * l(type_space[i]) / den
    end

    # new_u = 1 .- (sum(M, dims = 1) ./ n_divisions)[:]

    return new_u
end


function iterate_u(old_u::Vector, M::Matrix, tol::Float64 = 1e-10)
    res = Inf; # Initialize the residual to Inf
    u_last = old_u; # Save the last unmatched rates
    while res > tol
        # Calcualte the inegral in the denominator
        u_next = get_new_u(u_last, M)
        # calculate the residual
        res = norm(u_next[:] .- u_last[:])
        u_last = copy( u_next )
    end
    # u_last = 1 .- (sum(M, dims = 1) ./ n_divisions)

    return u_last[:]
end


function new_w(old_w::Vector, u::Vector, M::Matrix)
    # Calcualte the inegral in the denominator
    new_w = zeros(n_divisions)
    for i ∈ 1:n_divisions
        den = 1 + (θ * sum( M[i,:] .* u))
        num = θ * sum( M[i,:] .* (F[i, :] .- old_w) .* u)
        new_w[i] = num / den
    end
    return new_w
end


function iterate_w(old_w::Vector, u::Vector, M::Matrix; tol::Float64 = 1e-10)
    res = Inf; # Initialize the residual to Inf
    w_last = old_w; # Save the last unmatched rates
    while res > tol
        # Calcualte the inegral in the denominator
        w_next = new_w(w_last, u, M)
        # calculate the residual
        res = norm(w_next[:] .- w_last[:])
        # @show res
        w_last = copy( w_next )
    end
    return w_last[:]
end

function new_M(old_M::Matrix, u::Vector, w::Vector)
    
    S = zeros(n_divisions, n_divisions)
    for i ∈ 1:n_divisions, j ∈ 1:n_divisions
        S[i,j] = F[i,j] - w[i] - w[j]
    end
    
    return Int.(S .≥ 0)
    
end

function iterate_M(initial_M::Matrix, initial_u::Vector, initial_w::Vector, tol::Float64 = 1e-10; pLot::Bool = false)
    # Starting with initial guess of matching set, unmatched rates and value function iterate
    res = Inf; # Initialize the residual to Inf
    M_last = initial_M; # Save the last matching set
    u_last = initial_u; # Save the last unmatched rates
    w_last = initial_w; # Save the last value function
    if pLot
    println(
        UnicodePlots.heatmap(M_last, xfact=.1, yfact=.1, xoffset=-1.5,  legend = false,
        title = "Matching Set residual = $(res) iter = 0")#, colormap=:inferno)
        )
    end
    i = 1; # Initialize the iteration counter
    while res > tol
        
        # Calcualte the next unmatched rates
        u_next = iterate_u(u_last, M_last)
        # Calculate the next value function
        w_next = iterate_w(w_last, u_next, M_last, tol = 1e-4)
        # Calculate the next matching set
        M_next = new_M(M_last, u_next, w_next)
        # calculate the residual
        res = norm(M_next[:] .- M_last[:])
        if pLot
        println(
            UnicodePlots.heatmap(M_next, xfact=.1, yfact=.1, xoffset=-1.5, legend = false,
                    title = "Matching Set residual = $(res) iter = $i")#, colormap=:inferno)
                )
        else 
            println("Iteration $(i) residual = $(res)")
        end
        # Update the last values
        M_last = copy( M_next )
        u_last = copy( u_next )
        w_last = copy( w_next )
        i += 1
    end
    return M_last
end
