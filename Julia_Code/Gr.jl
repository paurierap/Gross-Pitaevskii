function Gr(M::SparseMatrixCSC{Float64,Int64},S::SparseMatrixCSC{Float64,Int64},Time_Simplices::Array{Float64,2},u0::Array{ComplexF64,1},Nodes::Array{Float64,1},Simplices::Array{Int64,2},Mesh2Space::Array{Int64,1},SpaceSize::Int64,k::Int64,r::Int64,β::Float64)

    # G_r(k) method with k order of Legendre polynomial.

    if r < 0 || r > k
        error("r must be such that 0 <= r <= k.")
    elseif r == 0     # r = 0 is equivalent to dG(k).
        return dGalerkin(M,S,Time_Simplices,u0,Nodes,Simplices,Mesh2Space,SpaceSize,k,β)
    elseif r == 1     # r = 1 is equivalent to cGP(k).
        return cGP(M,S,Time_Simplices,u0,Nodes,Simplices,Mesh2Space,SpaceSize,k,β)
    end

    N = size(Time_Simplices,2)
    τ = Time_Simplices[2,1]-Time_Simplices[1,1] # timestep.

    u = u0
    coeffs = vcat(u,zeros(ComplexF64,k*SpaceSize))*sqrt(τ) # Vector with the coefficients of each Legendre polynomial.
    # At the very beginning, the iterative method needs an initial guess.

    D0 = floor(Int64,0.5*(r-1))
    D1 = floor(Int64,0.5*r)

    mass = zeros(N+1)
    energy = zeros(N+1)
    mass[1] = mass_integration(u,Nodes,Simplices)
    energy[1] = energy_integration(β,u,Nodes,Simplices)

    A = spzeros(ComplexF64,(k+1)*SpaceSize,(k+1)*SpaceSize) # Initialize sparse matrix A that solves Ax = b,
    # where x is the vector of Legendre coefficients that is found iteratevely.
    # A is mostly constant through iterations.

    Id = Matrix(I,SpaceSize,SpaceSize) # Identity matrix.

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

    # Now, continuity conditions at the beginning of each interval:

    for j = 0:k
        A[(k-r+1)*SpaceSize+1:((k-r+1)+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += sqrt((2*j+1)/τ)*(-1)^j*Id
    end

    # Finally, continuity conditions on the derivatives depending on D0 and D1:

    for q = 1:D0
        i = k-r+1+q
        side = "left" # Left of each interval
        for j = 0:k
            A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += LegDerivative(q,j,τ,side)*M*im - LegDerivative(q-1,j,τ,side)*S
        end
    end

    for q = 1:D1
        i = k-r+1+q+D0
        side = "right" # Right of each interval
        for j = 0:k
            A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] += LegDerivative(q,j,τ,side)*M*im - LegDerivative(q-1,j,τ,side)*S
        end
    end

    # Let the magic begin:

    for j = 1:N

        if D0 == D1
            b = vcat(zeros((k-r+1)*SpaceSize),u,zeros(2*D1*SpaceSize))
        else
            b = vcat(zeros((k-r+1)*SpaceSize),u,zeros((2*D1-1)*SpaceSize))
        end

        A_matrix = copy(A) # A needs to be copied otherwise its values are changed.

        coeffs = Fixed_Point(β,b,coeffs,A,Nodes,Simplices,Mesh2Space,SpaceSize,k,r,Time_Simplices[:,j]) # Solve for the coefficients in the iterative solver.

        A = A_matrix

        fill!(u,0)

        for i = 0:k
            u += sqrt((2*i+1)/τ)*coeffs[i*SpaceSize+1:(i+1)*SpaceSize] # Recover actual solution from coefficients.
        end

        mass[j+1] = mass_integration(u,Nodes,Simplices) # Find the mass at that time.
        energy[j+1] = energy_integration(β,u,Nodes,Simplices) # Find the energy at that time.
    end

    # Matrix with all continuous derivatives:

    u_tot = zeros(ComplexF64,length(u),D0+1)
    u_tot[:,1] = u;

    for υ = 1:D0
        for i = 0:k
            u_tot[:,υ+1] += LegDerivative(υ,i,τ,"right")*coeffs[i*SpaceSize+1:(i+1)*SpaceSize]
        end
    end

    return u_tot,mass,energy

end
