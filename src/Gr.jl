"""
    Gr(M, S, time_intervals, u0, nodes, elements, map, ndofs, k, r, β;
       tol=1e-9, maxiter=100)

Solve iM u′ = S u + βΓ(u)u with the general Galerkin method Gᵣ(k).
Returns the final field and its D0 continuous derivatives, mass, and energy.
Time intervals must be contiguous, increasing and uniformly spaced.
"""
function Gr(M::SparseMatrixCSC{Float64,Int}, S::SparseMatrixCSC{Float64,Int},
            time_intervals::Matrix{Float64}, u0::Vector{ComplexF64},
            nodes::Vector{Float64}, elements::Matrix{Int}, map::Vector{Int},
            ndofs::Int, k::Int, r::Int, β::Float64;
            tol::Float64=1e-10, maxiter::Int=100)
    @assert 1 ≤ r ≤ k "r must satisfy 0 < r ≤ k."
    @assert tol > 0 && isfinite(tol) && maxiter > 0
    @assert size(time_intervals,1) == 2 && size(time_intervals,2) > 0
    @assert size(M) == size(S) == (ndofs,ndofs)
    @assert length(u0) == ndofs && length(nodes) == ndofs+2
    @assert map == [0; collect(1:ndofs); 0] "Expected homogeneous endpoint Dirichlet DOFs."
    @assert all(isfinite, u0) && isfinite(β)
    N = size(time_intervals,2)
    τ = time_intervals[2,1]-time_intervals[1,1]
    @assert τ > 0
    time_tol = 100eps(Float64)*max(1.0, maximum(abs, time_intervals))
    @assert all(isapprox.(time_intervals[2,:]-time_intervals[1,:], τ; atol=time_tol, rtol=1e-12)) "Uniform time steps are required."
    @assert all(isapprox.(time_intervals[1,2:end], time_intervals[2,1:end-1]; atol=time_tol, rtol=0)) "Time intervals must be contiguous."

    D0 = (r-1)÷2
    u = copy(u0)
    coeffs = vcat(sqrt(τ).*u, zeros(ComplexF64, k*ndofs))
    mass, energy = zeros(N+1), zeros(N+1)
    mass[1] = mass_integration(β, u, nodes, elements)
    energy[1] = energy_integration(β, u, nodes, elements)

    # Assemble once; all basis functions are expressed on the same local interval.
    A = assemble_time_matrix(M, S, k, r, τ)
    cache = NonlinearCache(nodes, elements, map, ndofs, k, r, τ)
    linear_factor = β == 0 ? lu(A) : nothing
    b = zeros(ComplexF64, (k+1)*ndofs)
    continuity = block_inds(k-r+1, ndofs)
    endpoint = [legendre_derivative_endpoint(0, s, τ) for s in 0:k]

    # Let the magic begin:
    for step in 1:N
        b[continuity] .= u
        if β == 0
            coeffs = linear_factor \ b
        else
            coeffs = fixed_point_iteration(β, b, coeffs, A, cache; tol, maxiter, step)
        end
        mul!(u, reshape(coeffs, ndofs, k+1), endpoint)
        mass[step+1] = mass_integration(β, u, nodes, elements)
        energy[step+1] = energy_integration(β, u, nodes, elements)
    end

    # Matrix with all continuous derivatives at the final right endpoint:
    E = [legendre_derivative_endpoint(q, s, τ) for s in 0:k, q in 0:D0]
    return reshape(coeffs, ndofs, k+1)*E, mass, energy
end

@inline block_inds(block, ndofs) = (block*ndofs+1):((block+1)*ndofs)

"""
    assemble_time_matrix(M, S, k, r, τ)

Assemble the linear variational, value-continuity and endpoint collocation rows.
"""
function assemble_time_matrix(M, S, k, r, τ)
    ndofs = size(M,1)
    D0, D1 = (r-1)÷2, r÷2
    I, J, V = Int[], Int[], ComplexF64[]
    Im, Jm, Vm = findnz(M)
    Is, Js, Vs = findnz(S)
    for j in 0:k, i in 0:(k-r)
        if i == j
            append!(I, Is .+ i*ndofs); append!(J, Js .+ j*ndofs); append!(V, -Vs)
        elseif is_coupled(j,i)
            append!(I, Im .+ i*ndofs); append!(J, Jm .+ j*ndofs)
            append!(V, (2im/τ)*sqrt((2i+1)*(2j+1)).*Vm)
        end
    end
    for j in 0:k
        append!(I, (1:ndofs) .+ (k-r+1)*ndofs)
        append!(J, (1:ndofs) .+ j*ndofs)
        append!(V, fill(legendre_derivative_endpoint(0,j,τ; side=:left), ndofs))
    end
    for (side, count, offset) in ((:left,D0,0), (:right,D1,D0))
        for q in 1:count, j in 0:k
            i = k-r+1+offset+q
            append!(I, Im .+ i*ndofs); append!(J, Jm .+ j*ndofs)
            append!(V, im*legendre_derivative_endpoint(q,j,τ; side).*Vm)
            append!(I, Is .+ i*ndofs); append!(J, Js .+ j*ndofs)
            append!(V, -legendre_derivative_endpoint(q-1,j,τ; side).*Vs)
        end
    end
    return sparse(I, J, V, (k+1)*ndofs, (k+1)*ndofs)
end

"""
    fixed_point_iteration(β, b, coeffs, A0, cache; tol, maxiter, step)

Freeze the density and solve (A0 - βΓ) x = b until both the coefficient
update and the nonlinear backward error converge. Sparse direct solves avoid
an inner iterative tolerance overwhelming the requested nonlinear accuracy.
"""
function fixed_point_iteration(β, b, coeffs, A0, cache; tol=1e-9, maxiter=100, step=1)
    x = copy(coeffs)
    Γ = nonlinear_matrix(x, cache)
    for iteration in 1:maxiter
        x_new = (A0 - β*Γ) \ b
        all(isfinite, x_new) || error("Nonfinite iterate at time step $step, iteration $iteration.")
        Δx = norm(x_new-x)
        Γ = nonlinear_matrix(x_new, cache)
        if Δx ≤ tol * (1+norm(x_new))
            A = A0 - β*Γ
            residual = A*x_new-b
            # Row scaling accounts for the different units of derivative rows.
            scale = max.(abs.(A)*abs.(x_new) + abs.(b), 1.0)
            if maximum(abs.(residual)./scale) ≤ tol
                return x_new
            end
        end
        x = x_new
    end
    error("Fixed-point iteration did not converge at time step $step after $maxiter iterations; reduce the time step.")
end
