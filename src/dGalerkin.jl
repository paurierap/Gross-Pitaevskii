function dGalerkin(M,S,Time_Simplices::Array{Float64,2},u0::Array{Float64,1},Nodes,Simplices,Mesh2Space,SpaceSize,k::Int64,β)

    # dG(k) method with k order of Legendre polynomial.

    N = size(Time_Simplices)[2]
    τ = Time_Simplices[2,1]-Time_Simplices[1,1]

    u = zeros(ComplexF64,SpaceSize,N+1)
    mass = zeros(N+1)
    energy = zeros(N+1)

    u[:,1] = u0

    A = spzeros(ComplexF64,(k+1)*SpaceSize,(k+1)*SpaceSize)
    b = spzeros(ComplexF64,(k+1)*SpaceSize,1)

    for i = 0:k

        for j = 0:k

            A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += (-1)^(i+j)*sqrt((2*i+1)*(2*j+1))/τ*M*im

            if i == j

                A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += S

            elseif ismult(i,j-1)

                A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += 2*sqrt((2*i+1)*(2*j+1))/τ*M*im

            end

        end

    end

    for j = 1:N

        t_n1 = Time_Simplices[1,j]
        t_n = Time_Simplices[2,j]

        for i = 0:k

            b[i*SpaceSize+1:(i+1)*SpaceSize] = (-1)^i*sqrt((2*i+1)/τ)*M*im*u[:,j]

        end

        A_matrix = copy(A) # A needs to be copied otherwise its values are changed.

        coeffs = Fixed_Point(β,b,coeffs,A,Nodes,Simplices,Mesh2Space,SpaceSize,k,r,Time_Simplices[:,j]) # Solve for the coefficients in the iterative solver.

        A = A_matrix

        fill!(u,0)

        for i = 0:k

            u += sqrt((2*i+1)/τ)*coeffs[i*SpaceSize+1:(i+1)*SpaceSize] # Recover actual solution from coefficients.

        end

    end

    return u,mass,energy

end
