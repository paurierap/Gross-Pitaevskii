"""
    assemble_global_matrix(nodes, elements, mesh_to_space_map, ndofs, local_assembler!)

Generic assembly routine for 1D FEM global matrices using triplet (COO) accumulation.

The function `local_assembler!` must have the signature
`local_assembler!(Mloc, element_nodes, ξ, ω)` and fill `Mloc` in place.
"""
function assemble_global_matrix(nodes::Vector{Float64}, elements::Matrix{Int}, mesh_to_space_map::Vector{Int}, ndofs::Int, local_assembler!::Function)

    nloc, nelem = size(elements)

    # Precompute reference quadrature points and weights
    ξ, ω = legendre_gauss_lobatto(nloc)

    # Preallocate storage. I, J and V form triplets (row,col,value) for a sparse matrix.
    I = Int[]
    J = Int[]
    V = Float64[]
    Mloc = zeros(nloc, nloc)

    for el = 1:nelem
        element = elements[:, el]
        element_nodes = nodes[element]

        # Fill local matrix in place
        local_assembler!(Mloc, element_nodes, ξ, ω)

        global_idx = mesh_to_space_map[element]

        for i = 1:nloc
            for j = 1:nloc
                ig, jg = global_idx[i], global_idx[j]
                if ig != 0 && jg != 0
                    push!(I, ig)
                    push!(J, jg)
                    push!(V, Mloc[i, j])
                end
            end
        end
    end

    return sparse(I, J, V, ndofs, ndofs)
end

"""
    assemble_mass_matrix(nodes, elements, mesh_to_space_map, ndofs)

Assemble the global mass matrix for a 1D finite element mesh.

# Arguments
- `nodes::Vector{Float64}` — global node coordinates.
- `elements::Matrix{Int}` — element connectivity, where `elements[:, el]`
  lists the node indices of element `el`.
- `mesh_to_space_map::Vector{Int}` — maps local node indices to global DOF indices.
- `ndofs::Int` — total number of degrees of freedom in the global system.

# Returns
- `M::SparseMatrixCSC{Float64, Int}` — assembled global mass matrix.
"""
assemble_mass_matrix(nodes, elements, map, ndofs) = assemble_global_matrix(nodes, elements, map, ndofs, assemble_mass_matrix_element!)

"""
    assemble_stiffness_matrix(nodes, elements, mesh_to_space_map, ndofs)

Assemble the global stiffness matrix for a 1D finite element mesh.

# Arguments
- `nodes::Vector{Float64}` — global node coordinates.
- `elements::Matrix{Int}` — element connectivity, where `elements[:, el]`
  lists the node indices of element `el`.
- `mesh_to_space_map::Vector{Int}` — maps local node indices to global DOF indices.
- `ndofs::Int` — total number of degrees of freedom in the global system.

# Returns
- `S::SparseMatrixCSC{Float64, Int}` — assembled global stiffness matrix.
"""
assemble_stiffness_matrix(nodes, elements, map, ndofs) = assemble_global_matrix(nodes, elements, map, ndofs, assemble_stiffness_matrix_element!)

"""
    assemble_mass_matrix_element!(Mloc::Matrix{Float64}, global_nodes::Vector{Float64}, ξ::Vector{Float64}, ω::Vector{Float64})

Compute the local mass matrix on a single 1D element whose physical
node coordinates are in `global_nodes`. Integration is performed using Legendre–Gauss–Lobatto quadrature.

# Arguments
- `Mloc::Matrix{Float64}`: preallocated local matrix.
- `global_nodes::Vector{Float64}`: global node coordinates.
- `ξ::Vector{Float64}`: Gauss-Lobato quadrature nodes in [-1,1].
- `ω::Vector{Float64}`: Gauss-Lobato quadrature weights in [-1,1].

# Returns
- `Mloc::Matrix{Float64}` — assembled local mass matrix.
"""
function assemble_mass_matrix_element!(Mloc::Matrix{Float64}, global_nodes::Vector{Float64}, ξ::Vector{Float64}, ω::Vector{Float64})

    nloc = length(global_nodes)
    h = abs(global_nodes[end] - global_nodes[1])

    # Mapping from local (ξ) to global coordinates
    xg = 0.5 * h * ξ .+ 0.5 * (global_nodes[end] + global_nodes[1])

    # Lagrange polynomials evaluated at xg
    φ = lagrange_matrix(global_nodes, xg)

    fill!(Mloc, 0.0)
    for i = 1:nloc
        for j = i:nloc  
            val = 0.5 * h * sum(ω .* φ[:, i] .* φ[:, j])
            Mloc[i, j] = val
            Mloc[j, i] = val
        end
    end

    return Mloc
end

"""
    assemble_stiffness_matrix_element!(Sloc::Matrix{Float64}, global_nodes::Vector{Float64}, ξ::Vector{Float64}, ω::Vector{Float64})

Compute the local stiffness matrix on a single 1D element whose physical
node coordinates are in `global_nodes`. Integration is performed using Legendre–Gauss–Lobatto quadrature.

# Arguments
- `Mloc::Matrix{Float64}`: preallocated local matrix.
- `global_nodes::Vector{Float64}`: global node coordinates.
- `ξ::Vector{Float64}`: Gauss-Lobato quadrature nodes in [-1,1].
- `ω::Vector{Float64}`: Gauss-Lobato quadrature weights in [-1,1].

# Returns
- `Sloc::Matrix{Float64}` — assembled local stiffness matrix.
"""
function assemble_stiffness_matrix_element!(Sloc::Matrix{Float64}, global_nodes::Vector{Float64}, ξ::Vector{Float64}, ω::Vector{Float64})
    nloc = length(global_nodes)
    h = abs(global_nodes[end] - global_nodes[1])
    
    # Mapping from local (ξ) to global coordinates
    xg = 0.5 * h * ξ .+ 0.5 * (global_nodes[end] + global_nodes[1])
    
    # Derivative of Lagrange polynomials evaluated at xg
    dφ = lagrange_derivative_matrix(global_nodes, xg)

    fill!(Sloc, 0.0)
    for i = 1:nloc
        for j = i:nloc
            val = 0.5 * h * sum(ω .* dφ[:, i] .* dφ[:, j])
            Sloc[i, j] = val
            Sloc[j, i] = val
        end
    end

    return Sloc
end