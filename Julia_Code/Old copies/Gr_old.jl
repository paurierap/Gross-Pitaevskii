function Gr(M,S,f,Time_Simplices::Array{Float64,2},u0::Array{Float64,1},Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,k::Int64,r::Int64,β::Float64)

    # G_r(k) method with k order of Legendre polynomial.

    if r < 0 || r > k

        error("r must be such that 0 <= r <= k.")

    elseif r == 0       # r = 0 is equivalent to dG(k).

        return dGalerkin(M,S,f,Time_Simplices,u0,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,k,β)

    elseif r == 1       # r = 1 is equivalent to cGP(k).

        return cGP(M,S,f,Time_Simplices,u0,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,k,β)

    end

    N = size(Time_Simplices,2)
    τ = Time_Simplices[2,1]-Time_Simplices[1,1] # timestep.

    u = u0
    coeffs = vcat(u0,zeros(k*SpaceSize))*sqrt(τ) # Vector with the coefficients of each Legendre polynomial.
    # At the very beginning, the iterative method needs an initial guess.

    D0 = floor(Int,0.5*(r-1))
    D1 = floor(Int,0.5*r)

    mass = zeros(N+1)
    energy = zeros(N+1)
    mass[1] = mass_integration(u,Nodes,Simplices,N_of_Simplices)
    energy[1] = energy_integration(β,u,Nodes,Simplices,N_of_Simplices)

    A = spzeros(ComplexF64,(k+1)*SpaceSize,(k+1)*SpaceSize) # Initialize sparse matrix A that solves Ax = b,
    # where x is the vector of Legendre coefficients that is found iteratevely.
    # A is mostly constant through iterations.

    Id = Matrix(I,SpaceSize,SpaceSize) # Identity matrix.

    # Auxiliar function to compute the scalar product P'j*Pi:

    function ismult(i,j)
        if j < 0 | i > j
            return false
        elseif i == j
            return true
        else
            ismult(i,j-2)
        end
    end

    # Double factorial function:

    function doublefact(x)
        y = 1;
        for j = 1:2:(2*x-1)
            y *= j
        end
        return y
    end

    # Construct matrix A. Let's begin with the local problem:

    for j = 0:k

        for i = 0:(k-r)

            if i == j

                A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] -= S # orthonormality.

            elseif ismult(i,j-1)

                A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += 2*sqrt((2*i+1)*(2*j+1))/τ*M*im

            end

        end

    end

    # Now continuity conditions at the beginning of each interval:

    for j = 0:k

        A[(k-r+1)*SpaceSize+1:((k-r+1)+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += sqrt((2*j+1)/τ)*(-1)^(j)*Id

    end

    # Finally, continuity conditions on the derivatives depending on D0 and D1:

    for q = 1:D0

        i = k-r+1+q

        A[i*SpaceSize+1:(i+1)*SpaceSize,1:SpaceSize] -= S/sqrt(τ)

        for j = 1:k

            A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += 1/τ*sqrt(2*j+1)* ( 2*(-1)^(j+q)*doublefact(q)*binomial(j+q,j-q)*M*im + (-1)^(j+q-1)*doublefact(q-1)*binomial(j+q-1,j-q-1)*S )

        end

    end

    for q = 1:D1

        i = k-r+1+q+D0

        A[i*SpaceSize+1:(i+1)*SpaceSize,1:SpaceSize] -= S/sqrt(τ);

        for j = 1:k

            A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += 1/τ*sqrt(2*j+1)* ( 2*doublefact(q)*binomial(j+q,j-q)*M*im + doublefact(q-1)*binomial(j+q-1,j-q-1)*S )

        end

    end

    for j = 1:N

        t_n1 = Time_Simplices[1,j]
        t_n = Time_Simplices[2,j]

        F = RHS_Integration(t_n1,t_n,k-r,D0,D1,f,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,"cGP")
        # Time integration of the RHS.

        #f_j = Assemble_RHS(Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,f,t_n) # Adjust terms to couple with continuity equations.
        #f_j1 = Assemble_RHS(Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,f,t_n1)

        #b = vcat(F,u,f_j,f_j1)

        # The following is ONLY because f(x,t) = 0, but it has to be generalized according to the RHS on the continuity of the derivatives:

        if D0 == D1
            b = vcat(F,u,zeros(2*D1*SpaceSize))
        else
            b = vcat(F,u,zeros((2*D1-1)*SpaceSize))
        end

        A_matrix = copy(A) # A needs to be copied otherwise its values are changed.

        coeffs = Fixed_Point(β,b,coeffs,A,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,k,r,Time_Simplices[:,j]) # Solve for the coefficients in the iterative solver.

        A = A_matrix

        fill!(u,0)

        for i = 0:k

            u += sqrt((2*i+1)/τ)*coeffs[i*SpaceSize+1:(i+1)*SpaceSize] # Recover actual solution from coefficients.

        end

        mass[j+1] = mass_integration(u,Nodes,Simplices,N_of_Simplices) # Find the mass at that time.
        energy[j+1] = energy_integration(β,u,Nodes,Simplices,N_of_Simplices) # Find the energy at that time.

    end

    return u,mass,energy

end
