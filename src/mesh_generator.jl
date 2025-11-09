"""
    generate_1D_mesh(x_left::Float64, x_right::Float64, Nx::Int, d::Int=1, t_start::Float64 = 0, t_end::Float64, N_t::Int, spacing::Symbol=:uniform)

Generate 1D mesh given the domain endpoints and a total number of nodes, and time discretization.

# Arguments
- `x_left::Float64`: left endpoint of the 1D domain.
- `x_right::Float64`: right endpoint of the 1D domain.
- `t_start::Float64 = 0`: starting time for a simulation.
- `t_end::Float64`: end time for a simulation.
- `d::Int64`: order of the Lagrange polynomials to be used as basis functions.
- `N_x::Int64`: number of space simplices in the mesh (equal to the number of nodes for `d=1`). 
- `N_t::Int64`: number of time steps.
- `spacing::Symbol=:uniform`: uniform spacing of the nodes by default Inverted chebyshev nodes (denser around the center of the domain) may be used with `:chebyshev`.

# Returns
`Nodes::Vector{Float64}` containing all nodes, `Simplices::Matrix{Float64}` containing nodes and element connect

# Notes
The order of the FEM discretization `d` directly influences the number of nodes, given that Lagrange polynomials are used as basis functions. Thus, the number of nodes increases for each simplex. Knowing the shape of the solution, we choose a spatial mesh with more nodes at the center, i.e. around 0, by means of traslation of Chebyshev nodes.
"""
function generate_mesh_1D(x_left::Float64, x_right::Float64, Nx::Int, t_start::Float64, t_end::Float64, Nt::Int; d::Int=1, spacing::Symbol=:uniform)

    N_nodes = Nx * d + 1

    if spacing == :chebyshev
        nodes = chebyshev_center_nodes(x_left, x_right, N_nodes)
    else
        nodes = collect(LinRange(x_left, x_right, N_nodes))
    end

    elements = zeros(Int, d+1, Nx) 
    for element = 1:Nx  
        elements[:,element] .= d * (element-1) .+ (1:(d+1))
    end

    # Time intervals
    t = collect(LinRange(t_start, t_end, Nt+1))
    time_intervals = Matrix([t[1:end-1] t[2:end]]')

    return nodes, elements, t, time_intervals
end

"""
    chebyshev_center_nodes(x_left::Float64, x_right::Float64, N::Int)

Generate `N` Chebyshev-type nodes in `[x_left, x_right]` that are denser
toward the **center** of the interval (as opposed to standard Chebyshev
nodes, which are denser near the boundaries).

The algorithm builds two mirrored Chebyshev half-domains and keeps
only the nodes within `[x_left, x_right]`.
"""
function chebyshev_center_nodes(x_left::Float64, x_right::Float64, N::Int)
    aL = (3 * x_left - x_right) * 0.5
    bL = (x_left + x_right) * 0.5
    aR = bL
    bR = (3 * x_right - x_left) * 0.5

    nodes_left = 0.5 .* ((bL - aL) .* cos.(π .* (N/2:-1:0) ./ N) .+ aL .+ bL)
    nodes_right = 0.5 .* ((bR - aR) .* cos.(π .* (N-1:-1:N/2) ./ N) .+ aR .+ bR)

    nodes = vcat(nodes_left, nodes_right)
    nodes = filter(x -> x ≥ x_left && x ≤ x_right, nodes)
    
    return nodes
end