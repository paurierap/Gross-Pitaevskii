#= Gross-Pitaevskii equation solver for general beta:

id_t(u)=-d_xx(u)+β|u|^2u

=#

using LinearAlgebra, SparseArrays, JLD, Plots, Debugger

include("LoadCode.jl")
@time begin
left_boundary = -20.
right_boundary = 20.
β = -2.

Nx = 2^10
Nt = 2^10

bw = load("1dMesh_"*string(Nx)*"x"*string(Nt)*".jld")
Nodes = bw["Nodes"]
Simplices=convert(Array{Int64,2}, bw["Simplices"])
Time_Nodes = bw["Time_Nodes"]
Time_Simplices = bw["Time_Simplices"]
N_of_Nodes = length(Nodes)
N_of_Simplices = size(Simplices,2)

Mesh2Space, SpaceSize = Mixed_BC(left_boundary,right_boundary,Nodes,N_of_Nodes,"Dirichlet")

# Assemble stiffness matrix:

S = Assemble_S(Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices)

# Assemble mass matrix:

M = Assemble_M(Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices)

# General Galerkin:

f(t,x) = 0

x = Nodes[2:end-1] # Mesh without first and last nodes.

u0 = (8*(9*exp.(-4*x)+16*exp.(4*x))-32*(4*exp.(-2*x)+9*exp.(2*x)))./(-128 .+4*exp.(-6*x)+16*exp.(6*x)+81*exp.(-2*x)+64*exp.(2*x)) # Initial condition.

k = 8 # G_r(k).
r = 4
end
error()
u,mass,energy = Gr(M,S,f,Time_Simplices,u0,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,k,r,β)

# Exact solution u(x,1):

sol(x,t) = (8*exp(4im*t)*(9*exp.(-4*x)+16*exp.(4*x))-32*exp(16im*t)*(4*exp.(-2*x)+9*exp.(2*x)))./(-128*cos(12*t) .+4*exp.(-6*x)+16*exp.(6*x)+81*exp.(-2*x)+64*exp.(2*x))

#p = plot(x,abs.(u))
#plot!(x,abs.(sol),linestyle = :dash)
