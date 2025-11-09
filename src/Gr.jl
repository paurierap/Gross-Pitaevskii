function Gr(M::SparseMatrixCSC{Float64,Int64}, 
            S::SparseMatrixCSC{Float64,Int64},
            time_intervals::Matrix{Float64}, 
            u0::Vector{ComplexF64}, 
            nodes::Vector{Float64}, 
            elements::Matrix{Int}, 
            mesh_to_space_map::Vector{Int}, 
            ndofs::Int64, k::Int64, r::Int64, β::Float64)

    @assert 1 ≤ r ≤ k "r must satisfy 0 < r ≤ k."

    N_time_intervals = size(time_intervals, 2)
    τ = time_intervals[2,1] - time_intervals[1,1] # timestep

    u = u0
    coeffs = vcat(u, zeros(ComplexF64, k * ndofs)) * sqrt(τ) 

    D0 = floor(Int64,0.5*(r-1))
    D1 = floor(Int64,0.5*r)

    mass = zeros(N_time_intervals+1)
    energy = zeros(N_time_intervals+1)
    mass[1] = mass_integration(β, u, nodes, elements)
    energy[1] = energy_integration(β, u, nodes, elements)

    # Build A using triplet assembly:
    I, J, V = Int[], Int[], ComplexF64[]
    Im, Jm, Vm = findnz(M)
    Is, Js, Vs = findnz(S)

    # Mapping for block indices:
    #block_inds(i, block) = (i*block + 1):((i+1)*block)

    #= Construct matrix A. Let's begin with the local problem:
    for j = 0:k
        col_block = block_inds(j, ndofs)
        for i = 0:(k-r)
            row_block = block_inds(i, ndofs)
            if i == j
                append!(I, repeat(row_block, outer=ndofs))
                append!(J, repeat(col_block, inner=ndofs))
                append!(V, vec(-S))
            elseif is_coupled(i,j)
                append!(I, repeat(row_block, outer=ndofs))
                append!(J, repeat(col_block, inner=ndofs))
                append!(V, vec((2/τ) * sqrt((2*i+1) * (2*j+1)) * im * M))
            end
        end
    end

    # Now, continuity conditions at the beginning of each interval:
    row_block = (k - r + 1) * ndofs
    for j = 0:k
        col_block = j * ndofs
        append!(I, row_block .+ (1:ndofs))
        append!(J, col_block .+ (1:ndofs))
        append!(V, fill((-1)^j * sqrt((2*j + 1) / τ), ndofs))
    end

    #Finally, continuity conditions on the derivatives depending on D0 and D1:
    for q = 1:D0
        row_block = block_inds(k - r + 1 + q, ndofs)
        for j = 0:k
            col_block = block_inds(j, ndofs)
            append!(I, repeat(row_block, outer=ndofs))
            append!(J, repeat(col_block, inner=ndofs))
            append!(V, vec(legendre_derivative_endpoint(q, j, τ, side=:left) * im * M - legendre_derivative_endpoint(q-1, j, τ, side=:left) * S))
        end
    end

    for q = 1:D1
        row_block = block_inds(k - r + 1 + q + D0, ndofs)
        for j = 0:k
            col_block = block_inds(j, ndofs)
            append!(I, repeat(row_block, outer=ndofs))
            append!(J, repeat(col_block, inner=ndofs))
            append!(V, vec(legendre_derivative_endpoint(q, j, τ, side=:right) * im * M - legendre_derivative_endpoint(q-1, j, τ, side=:right) * S))
        end
    end
    =#

    for j = 0:k
        for i = 0:(k-r)
            if i == j
                append!(I, Is .+ i * ndofs)
                append!(J, Js .+ j * ndofs)
                append!(V, -Vs)
            elseif is_coupled(j,i)
                append!(I, Im .+ i * ndofs)
                append!(J, Jm .+ j * ndofs)
                append!(V, (2/τ) * sqrt((2*i+1) * (2*j+1)) * im * Vm)
            end
        end
    end
    # Now, continuity conditions at the beginning of each interval:
    i = k - r + 1
    for j = 0:k
        append!(I, collect(1:ndofs) .+ i * ndofs)
        append!(J, collect(1:ndofs) .+ j * ndofs)
        append!(V, fill((-1)^j * sqrt((2*j + 1) / τ), ndofs))
    end

    #Finally, continuity conditions on the derivatives depending on D0 and D1:
    
    for q = 1:D0
        i = k - r + 1 + q
        for j = 0:k
            append!(I, Im .+ i * ndofs)
            append!(J, Jm .+ j * ndofs)
            append!(V, Vm * legendre_derivative_endpoint(q, j, τ, side=:left) * im)
            append!(I, Is .+ i * ndofs)
            append!(J, Js .+ j * ndofs)
            append!(V, -Vs * legendre_derivative_endpoint(q-1, j, τ, side=:left))
        end
    end

    
    for q = 1:D1
        i = k - r + 1 + q + D0
        for j = 0:k
            append!(I, Im .+ i * ndofs)
            append!(J, Jm .+ j * ndofs)
            append!(V, Vm * legendre_derivative_endpoint(q, j, τ, side=:right) * im)
            append!(I, Is .+ i * ndofs)
            append!(J, Js .+ j * ndofs)
            append!(V, -Vs * legendre_derivative_endpoint(q-1, j, τ, side=:right))
        end
    end

    A = sparse(I, J, V, (k+1) * ndofs, (k+1) * ndofs)
    
    # Let the magic begin:
    for j = 1:N_time_intervals 
        if D0 == D1
            b = vcat(zeros((k-r+1) * ndofs), u, zeros(2 * D1 * ndofs))
        else
            b = vcat(zeros((k-r+1) * ndofs), u, zeros((2 * D1 - 1) * ndofs))
        end

        # Use LU factorization of A as preconditioner for iterative solvers:
        A_fact = lu(A) 

        # Solve for the coefficients in the iterative solver:
        coeffs = fixed_point_iteration(β, b, coeffs, A, A_fact, nodes, elements,mesh_to_space_map, ndofs, k, r, time_intervals[:,j]) 

        fill!(u,0)
        for i = 0:k
            u += sqrt((2*i+1) / τ) * coeffs[block_inds(i, ndofs)] 
        end

        mass[j+1] = mass_integration(β, u, nodes, elements)
        energy[j+1] = energy_integration(β, u, nodes, elements)
    end

    # Matrix with all continuous derivatives:
    u_tot = zeros(ComplexF64,length(u),D0+1)
    u_tot[:,1] = u;

    for υ = 1:D0
        for i = 0:k
            u_tot[:,υ+1] += LegDerivative(υ,i,τ,"right")*coeffs[i*ndofs+1:(i+1)*ndofs]
        end
    end

    return u_tot,mass,energy

end

# Iterative algorithm to calculate the solution for a time interval. Usually converges in 6-8 iterations!

function fixed_point_iteration(β::Float64, 
                     b::Vector{ComplexF64},
                     x_prev::Vector{ComplexF64},
                     A0::SparseMatrixCSC{ComplexF64,Int64},
                     A0_fact::SparseArrays.UMFPACK.UmfpackLU{ComplexF64},
                     nodes::Vector{Float64},
                     elements::Matrix{Int64},
                     mesh_to_space_map::Vector{Int64},
                     ndofs::Int64,
                     k::Int64, r::Int64, 
                     time_interval::Vector{Float64};
                     tol::Float64 = 1e-9, 
                     tol_gmres::Float64 = 1e-6, 
                     maxiter::Int64 = 100)

    τ = abs(time_interval[2]-time_interval[1]) # timestep
    
    # Preallocate triplet containers (grow as needed but reuse)
    I = Int[]
    J = Int[]
    V = ComplexF64[]

    # Collocation method parameters
    D0 = floor(Int,0.5*(r-1))
    D1 = floor(Int,0.5*r)
    
    Δx = 1.0
    it = 0

    while Δx > tol
        it += 1

        if it == maxiter
            error("Maximum number of iterations reached. Solution has not converged")
        end

        # Build ΔA(x_prev) in triplet form
        empty!(I); empty!(J); empty!(V)

        # Loop over block indices (use same indexing as before)
        for i = 0:(k - r)
            for j = 0:k
                Ig, Jg, Vg = integrate_Γ(x_prev, i, j, k, time_interval, nodes, elements, mesh_to_space_map, ndofs)
                append!(I, Ig .+ i * ndofs)
                append!(J, Jg .+ j * ndofs)
                append!(V, -β * Vg)
            end
        end
        
        for q = 1:D0
            i = k - r + 1 + q
            for j = 0:k
                for α = 0:(q - 1)
                    Ig, Jg, Vg = derivative_Γ(q-1-α, k, τ, x_prev, nodes, elements, mesh_to_space_map, ndofs, side=:left)
                    append!(I, Ig .+ i * ndofs)
                    append!(J, Jg .+ j * ndofs)
                    append!(V, -β * Vg * binomial(q-1,α) * legendre_derivative_endpoint(α, j, τ, side=:left))
                end
            end
        end
        
        for q = 1:D1
            i = k - r + 1 + q + D0
            for j = 0:k
                for α = 0:(q - 1)
                    Ig, Jg, Vg = derivative_Γ(q-1-α, k, τ, x_prev, nodes, elements, mesh_to_space_map, ndofs, side=:left)
                    append!(I, Ig .+ i * ndofs)
                    append!(J, Jg .+ j * ndofs)
                    append!(V, -β * Vg * binomial(q-1,α) * legendre_derivative_endpoint(α, j, τ, side=:right))
                end
            end
        end 
        
        # Build ΔA as sparse
        ΔA = sparse(I, J, V, (k+1) * ndofs, (k+1) * ndofs)

        A = A0 + ΔA
        x_new = gmres(A, b; abstol = tol_gmres) 

        if norm(A * x_new - b) / norm(b) > 1e-5
            x_new = A \ b
        end

        Δx = norm(x_new - x_prev)
        x_prev = x_new

        @info "FixedPoint iteration $it, Δx = $Δx"
    end

    return x_new
end