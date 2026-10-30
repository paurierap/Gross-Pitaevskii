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

    @assert x_left < x_right "Spatial endpoints must be increasing."
    @assert t_start < t_end "Time endpoints must be increasing."
    @assert Nx > 0 && Nt > 0 && d > 0 "Mesh counts and FEM degree must be positive."
    @assert spacing in (:uniform, :chebyshev) "Unknown mesh spacing."

    # Subdivide each physical element into d equal parts.
    vertices = spacing == :chebyshev ? chebyshev_center_nodes(x_left, x_right, Nx + 1) :
               collect(LinRange(x_left, x_right, Nx + 1))
    nodes = Vector{Float64}(undef, Nx * d + 1)
    for el in 1:Nx, j in 0:(d-1)
        nodes[(el-1)*d+j+1] = vertices[el] + (j/d)*(vertices[el+1]-vertices[el])
    end
    nodes[end] = x_right

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

Generate `N` cosine-spaced nodes in `[x_left, x_right]` that are denser
toward the **center** of the interval (as opposed to standard Chebyshev
nodes, which are denser near the boundaries).

The mapping `±(1 - cos(π|s|/2))`, with `s` uniformly spaced in `[-1, 1]`,
preserves both endpoints and is symmetric about the midpoint.
"""
function chebyshev_center_nodes(x_left::Float64, x_right::Float64, N::Int)
    @assert N ≥ 2 && x_left < x_right
    midpoint = (x_left + x_right) / 2
    radius = (x_right - x_left) / 2
    s = range(-1.0, 1.0; length=N)
    nodes = [midpoint + sign(z)*radius*(1-cospi(z/2)) for z in s]
    nodes[1], nodes[end] = x_left, x_right
    return nodes
end
