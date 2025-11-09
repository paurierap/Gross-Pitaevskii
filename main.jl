"""
    Gross-Pitaevskii equation solver for general β:
			i ∂u/∂t = -∂²u/∂x² + β|u|²u
		
Here, this is specialized for the two interacting solitons problem in 1D, with:
- β = -2,
- x ∈ R,
- t ∈ [0,T), 
- Dirichlet boundary conditions at ±∞, ie u(-∞) = u(+∞) = 0.
- Initial condition: 
u(0, x) = (8 * (9*exp(-4*x) + 16*exp(4*x)) - 32*(4*exp(-2*x) + 9*exp(2*x)))/(-128 + 4*exp(-6*x) + 16*exp(6*x) +81*exp(-2*x) +64*exp.(2*x)). 

Provided that u(0, x) shows asymptotic behaviour at ±∞ towards 0, the BC's can be enforced much closer to the origin.

This problem has the analytical solution: 
u(x, t) = (8*exp(4it) * (9*exp(-4*x) + 16*exp(4*x)) - 32*exp(16it) * (4*exp(-2*x) + 9*exp(2*x)))/(-128 * cos(t) + 4*exp(-6*x) + 16*exp(6*x) +81*exp(-2*x) +64*exp.(2*x))

From here, we see that the mass and energy of the system are, respectively:
- M[u] = ∫ |u|² dx = 12,
- E[u] = ∫ |∂u/∂x|² + ½β|u|⁴ dx = -48.
"""

##########################################
##### Include dependencies and utils #####
##########################################

using LinearAlgebra, SparseArrays, JLD, IterativeSolvers
include("src/LoadCode.jl")

#############################
##### Define constants  #####
#############################

const left_boundary = -20.
const right_boundary = 20.
const t_start = 0.
const t_end = 5.
const d = 2 # Order of FEM basis functions.
const Nx = 1 << 8
const Nt = 1 << 9
const β = -2.
const k = 3 # G_r(k).
const r = 1

###################################
##### Spatial discretization  #####
###################################

# Compute mesh and time intervals:
nodes, elements, _ , time_intervals = generate_mesh_1D(left_boundary, right_boundary, Nx, t_start, t_end, Nt, spacing=:chebyshev)

# Apply Dirichlet boundary conditions:
ndofs = length(nodes) - 2
mesh_to_space_map = [0; collect(1:ndofs); 0]
x = nodes[2:end-1] 

# Assemble stiffness matrix:
S = assemble_stiffness_matrix(nodes, elements, mesh_to_space_map, ndofs)

# Assemble mass matrix:
M = assemble_mass_matrix(nodes, elements, mesh_to_space_map, ndofs)

#######################################
##### Solve with general Galerkin #####
#######################################

# Initial condition:
u0 = complex(float((8*(9*exp.(-4*x)+16*exp.(4*x))-32*(4*exp.(-2*x)+9*exp.(2*x)))./(-128 .+4*exp.(-6*x)+16*exp.(6*x)+81*exp.(-2*x)+64*exp.(2*x))))

# Start time integration
u, mass, energy = Gr(M, S, time_intervals, u0, nodes, elements, mesh_to_space_map,ndofs, k, r, β)

save("solution_file_"*string(Nx)*"x"*string(Nt)*"_k"*string(k)*"_r"*string(r)*"_P"*string(size(Simplices,1)-1)*"_T"*string(int(t_end))*".jld","u", u, "mass",mass, "energy", energy)