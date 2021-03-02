function cGP(M,S,Time_Simplices::Array{Float64,2},u0::Array{Float64,1},Nodes,Simplices,Mesh2Space,SpaceSize,k::Int64,β)

    # cGP(k) method with k order of Legendre polynomial.

    N = size(Time_Simplices)[2]
    τ = Time_Simplices[2,1]-Time_Simplices[1,1]

    #u = zeros(ComplexF64,SpaceSize,N+1);

    u = u0
    coeffs = vcat(u0,zeros(k*SpaceSize))*sqrt(τ)
    Id = Matrix(I,SpaceSize,SpaceSize)

    A = spzeros(ComplexF64,(k+1)*SpaceSize,(k+1)*SpaceSize)

    mass = zeros(N+1)
    energy = zeros(N+1)
    mass[1] = mass_integration(u,Nodes,Simplices)
    energy[1] = energy_integration(β,u,Nodes,Simplices)

    for i = 0:(k-1)

        for j = 0:k

            if i == j

                A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] -= S

            elseif ismult(i,j-1)

                A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += 2*sqrt((2*i+1)*(2*j+1))/τ*M*im

            end

        end

    end

    for j = 0:k

        A[k*SpaceSize+1:(k+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += (-1)^j*sqrt((2*j+1)/τ)*Id

    end

    for j = 1:N

        t_n1 = Time_Simplices[1,j]; t_n = Time_Simplices[2,j]

        b = vcat(zeros(k*SpaceSize),u)

        A_matrix = copy(A)

        coeffs = Fixed_Point(β,b,coeffs,A,Nodes,Simplices,Mesh2Space,SpaceSize,k,1,Time_Simplices[:,j])

        A = A_matrix

        fill!(u,0)

        for i = 0:k

            u += sqrt((2*i+1)/τ)*coeffs[i*SpaceSize+1:(i+1)*SpaceSize]

        end

        mass[j+1] = mass_integration(u,Nodes,Simplices)
        energy[j+1] = energy_integration(β,u,Nodes,Simplices)

    end

    return u,mass,energy

end
