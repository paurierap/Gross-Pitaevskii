#= Gross-Pitaevskii equation solver for general beta:

id_t(u)=-d_xx(u)+β|u|^2u

=#

using LinearAlgebra, SparseArrays, JLD

include("LoadCode.jl")

left_boundary = -20.
right_boundary = 20.
β = -2.

Nx = 2^6
Nt = 2^14
d = 1

bw = load("1dMesh_"*string(Nx)*"x"*string(Nt)*"_P"*string(d)*".jld")
Nodes = bw["Nodes"]
Simplices = convert(Array{Int64,2}, bw["Simplices"])
Time_Nodes = bw["Time_Nodes"]
Time_Simplices = bw["Time_Simplices"]

Mesh2Space, SpaceSize = Mixed_BC(left_boundary,right_boundary,Nodes,"Dirichlet")

# Assemble stiffness matrix:

S = Assemble_S(Nodes,Simplices,Mesh2Space,SpaceSize)

# Assemble mass matrix:

M = Assemble_M(Nodes,Simplices,Mesh2Space,SpaceSize)

# General Galerkin:

x = Nodes[2:end-1] # Mesh without first and last nodes.

u0 = (8*(9*exp.(-4*x)+16*exp.(4*x))-32*(4*exp.(-2*x)+9*exp.(2*x)))./(-128 .+4*exp.(-6*x)+16*exp.(6*x)+81*exp.(-2*x)+64*exp.(2*x)) # Initial condition.

k = 4 # G_r(k).
r = 3

u,mass,energy = Gr(M,S,Time_Simplices,u0,Nodes,Simplices,Mesh2Space,SpaceSize,k,r,β)

# Exact solution u(x,t):

# sol(x,t) = (8*exp(4im*t)*(9*exp.(-4*x)+16*exp.(4*x))-32*exp(16im*t)*(4*exp.(-2*x)+9*exp.(2*x)))./(-128*cos(12*t) .+4*exp.(-6*x)+16*exp.(6*x)+81*exp.(-2*x)+64*exp.(2*x))

save("sol_file_"*string(Nx)*"x"*string(Nt)*"_k"*string(k)*"_r"*string(r)*"_P"*string(size(Simplices,1)-1)*".jld","u",u, "mass",mass,"energy",energy)

#p = plot(x,abs.(u))
#plot!(x,abs.(sol),linestyle = :dash)
