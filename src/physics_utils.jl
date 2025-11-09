"""
    quantity_integration(β, u, nodes, elements, f)

Generic FEM quadrature-based integration routine for computing physical
quantities (e.g., energy, mass) in 1D systems.

The integrand function `f(β, u_loc, du_loc)` should return a *vector* of
values at quadrature points, representing the pointwise quantity to integrate.

# Arguments
- `β::Float64`: Nonlinearity parameter.
- `u::Vector{ComplexF64}`: Discrete solution vector (excluding Dirichlet BCs).
- `nodes::Vector{Float64}`: Node coordinates.
- `elements::Matrix{Int}`: Element connectivity.
- `f::Function`: Function defining the integrand at quadrature points.

# Returns
- `quantity::Float64`: Integrated scalar quantity.
"""
function quantity_integration(β::Float64, u::Vector{ComplexF64}, nodes::Vector{Float64}, elements::Matrix{Int}, f::Function)

    nloc, nelem = size(elements)
    quantity = 0.0

    # Extend solution for homogeneous Dirichlet BCs
    u = vcat(0.0 + 0im, u, 0.0 + 0im)

    # 5-point Gauss quadrature
    ξ = Quad_5.x
    ω = Quad_5.w

    # Preallocate matrices
    φ = zeros(length(ξ), nloc)
    dφ = similar(φ)

    for el in 1:nelem
        elem_nodes = nodes[elements[:, el]]
        h = abs(elem_nodes[end] - elem_nodes[1])

        u_el = u[elements[:, el]]

        # Map quadrature points to physical domain
        xq = 0.5 * h * ξ .+ 0.5 * (elem_nodes[end] + elem_nodes[1])

        # Evaluate basis and derivatives
        φ .= lagrange_matrix(elem_nodes, xq)
        dφ .= lagrange_derivative_matrix(elem_nodes, xq)

        # Field and derivative at quadrature points
        u_loc = φ * u_el
        du_loc = dφ * u_el

        # Integrate local contribution (f returns vector)
        quantity += 0.5 * h * sum(ω .* f(β, u_loc, du_loc))
    end

    return quantity
end

"""
    mass_integration(β::Float64, u::Vector{ComplexF64}, nodes::Vector{Float64}, elements::Matrix{Int})

Compute the total mass (or number of particles) of the 1D Gross–Pitaevskii system:
    M[u] = ∫ |u|² dx

Assumes Dirichlet boundary conditions u=0 at both ends.

# Arguments
- `β::Float64`: Gross-Pitaevskii interaction term.
- `u::Vector{ComplexF64}`: vectorized solution of the system.
- `nodes::Vector{Float64}`: global node coordinates.
- `elements::Matrix{Int}`: element connectivity, where `elements[:, el]`
  lists the node indices of element `el`.

# Returns
mass::{Float64}, the corresponding total mass.
"""
mass_integration(β, u, nodes, elements) =
    quantity_integration(β, u, nodes, elements, mass_integrand)

mass_integrand(β, u_q, du_q) = abs2.(u_q)

"""
    energy_integration(β::Float64, u::Vector{ComplexF64}, nodes::Vector{Float64}, elements::Matrix{Int})

Compute the total energy of the 1D Gross–Pitaevskii system:
    E[u] = ∫ |∂u/∂x|² + ½β|u|⁴ dx

Assumes Dirichlet boundary conditions u=0 at both ends.

# Arguments
- `β::Float64`: Gross-Pitaevskii interaction term.
- `u::Vector{ComplexF64}`: vectorized solution of the system.
- `nodes::Vector{Float64}`: global node coordinates.
- `elements::Matrix{Int}`: element connectivity, where `elements[:, el]`
  lists the node indices of element `el`.

# Returns
energy::{Float64}, the corresponding total energy.
"""
energy_integration(β, u, nodes, elements) =
    quantity_integration(β, u, nodes, elements, energy_integrand)

energy_integrand(β, u_q, du_q) = abs2.(du_q) .+ 0.5 * β .* abs2.(u_q).^2