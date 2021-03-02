function cGP(M,A,f,Time_Simplices::Array{Float64,2},u0::Array{Float64,1},Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,k::Int64,β)

    # cGP(k) method with k order of Legendre polynomial.

    N = size(Time_Simplices)[2]
    τ = Time_Simplices[2,1]-Time_Simplices[1,1]

    #u = zeros(ComplexF64,SpaceSize,N+1);

    u = u0
    coeffs = vcat(u0,zeros(k*SpaceSize))*sqrt(τ)
    Id = Matrix(I,SpaceSize,SpaceSize)

    B = spzeros(ComplexF64,(k+1)*SpaceSize,(k+1)*SpaceSize)

    mass = zeros(size(Time_Simplices)[2]+1)
    energy = zeros(size(Time_Simplices)[2]+1)
    mass[1] = mass_integration(u,Nodes,Simplices,N_of_Simplices)
    energy[1] = energy_integration(beta,u,Nodes,Simplices,N_of_Simplices)

    # Auxiliar function to compute the scalar product P'j*Pi.

    function ismult(i,j)
        if j < 0 | i > j
            return false
        elseif i == j
            return true
        else
            ismult(i,j-2)
        end
    end

    for i = 0:(k-1)

        for j = 0:k

            if i == j

                B[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] -= A

            elseif ismult(i,j-1)

                B[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += 2*sqrt((2*i+1)*(2*j+1))/τ*M*im

            end

        end

    end

    for j = 0:k

        B[k*SpaceSize+1:(k+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += (-1)^j*sqrt((2*j+1)/τ)*Id;

    end

    for j = 1:N

        t_n1 = Time_Simplices[1,j]; t_n = Time_Simplices[2,j];

        F = RHS_Integration(t_n1,t_n,k-1,0,0,f,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,"cGP")

        b = vcat(F,u);

        B_matrix = copy(B);

        coeffs = Fixed_Point(β,b,coeffs,B,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,k,1,Time_Simplices[:,j]);

        B = B_matrix;

        fill!(u,0)

        for i = 0:k

            u += sqrt((2*i+1)/τ)*coeffs[i*SpaceSize+1:(i+1)*SpaceSize]

        end

        mass[j+1] = mass_integration(u,Nodes,Simplices,N_of_Simplices)
        energy[j+1] = energy_integration(β,u,Nodes,Simplices,N_of_Simplices)

    end

    return u,mass,energy

end
