function dGalerkin(M,A,f,Time_Simplices::Array{Float64,2},u0::Array{Float64,1},Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,k::Int64,β)

    # dG(k) method with k order of Legendre polynomial.

    N = size(Time_Simplices)[2]
    τ = Time_Simplices[2,1]-Time_Simplices[1,1]

    u = zeros(ComplexF64,SpaceSize,N+1)
    mass = zeros(N+1)
    energy = zeros(N+1)

    u[:,1] = u0

    B = spzeros(ComplexF64,(k+1)*SpaceSize,(k+1)*SpaceSize)
    b = spzeros(ComplexF64,(k+1)*SpaceSize,1)

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

    for i = 0:k

        for j = 0:k

            B[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += (-1)^(i+j)*sqrt((2*i+1)*(2*j+1))/τ*M*im

            if i == j

                B[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += A

            elseif ismult(i,j-1)

                B[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += 2*sqrt((2*i+1)*(2*j+1))/τ*M*im

            end

        end

    end

    for j = 1:N

        t_n1 = Time_Simplices[1,j]
        t_n = Time_Simplices[2,j]

        for i = 0:k

            b[i*SpaceSize+1:(i+1)*SpaceSize] = (-1)^i*sqrt((2*i+1)/τ)*M*im*u[:,j]

        end

        F = RHS_Integration(t_n1,t_n,k,0,0,f,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,"dG")

        b += F

        u[:,j+1] = Fixed_Point(β,b,u[:,j],B,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,k,τ)

        mass[j] = mass_integration(u[:,j+1],Nodes,Simplices,N_of_Simplices)
        energy[j] = energy_integration(β,u[:,j+1],Nodes,Simplices,N_of_Simplices)

    end

    return u,mass,energy

end
