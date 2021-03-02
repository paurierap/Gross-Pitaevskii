# Generate 1D mesh and time discretization.

using JLD

x_left = -20.
x_right = 20.

t_start = 0
t_end = 30

d = 1 # Order of FEM basis functions.

N_x = 2^6
N_t = 2^19

Nodes_t = Array(LinRange(t_start,t_end,N_t+1)) # Time nodes equispaced over [0,1]

#Nodes = Array(LinRange(x_left,x_right,N_x+1))

# Knowing the shape of the solution, we choose a spatial mesh with more nodes at the center,
# i.e. around 0, by means of traslation of Chebyshev nodes.

a = (3*x_left-x_right)*0.5; b = (x_left+x_right)*0.5;
N = N_x*d
Nodes_left = 0.5*Array((b-a)*cos.(π*LinRange(N/2,0,convert(Int64,N/2+1))/N).+a.+b);
a = b; b = 0.5*(3*x_right-x_left);
Nodes_right = 0.5*Array((b-a)*cos.(π*LinRange(N-1,N/2,convert(Int64,N/2))/N).+a.+b);
Nodes = vcat(Nodes_left,Nodes_right);

# Space simplicial elements of P-FEM of order d:

Simplices = zeros(d+1,N_x)

for i = 1:N_x

    for j = 1:(d+1)

        Simplices[j,i] = d*(i-1)+j

    end

end

# Time simplicial elements:

Time_Simplices = zeros(2,N_t)

for i = 1:N_t

    Time_Simplices[1,i] = Nodes_t[i]
    Time_Simplices[2,i] = Nodes_t[i+1]

end

save("1dMesh_"*string(N_x)*"x"*string(N_t)*"_P"*string(d)*".jld", "Nodes", Nodes, "Simplices",Simplices,"Time_Nodes", Nodes_t, "Time_Simplices",Time_Simplices)
