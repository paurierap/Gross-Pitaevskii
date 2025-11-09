function integrate_Γ(coeffs::Vector{ComplexF64},
                      p::Int, q::Int, k::Int,
                      time_interval::Vector{Float64},
                      nodes::Vector{Float64},
                      elements::Matrix{Int},
                      mesh_to_space_map::Vector{Int},
                      ndofs::Int)

    nloc, nelem = size(elements)
    τ = abs(time_interval[2] - time_interval[1])

    # --- Quadrature setup ---

    # 5-point Gaussian quadrature in space
    ξ = Quad_5.x
    ω = Quad_5.w
    len_x = length(ξ)

    # (k+1)-point Legendre–Gauss–Lobatto in time
    η, ζ = legendre_gauss_lobatto(k + 1)
    t = 0.5 * τ .* η .+ 0.5 * (time_interval[1] + time_interval[2])

    # Legendre basis in time
    P = legendre_matrix(k, time_interval, t)

    # --- Reconstruct u(x, t) from Legendre coefficients ---
    # u has size (ndofs, len_t)
    len_t = length(t)
    u = zeros(ComplexF64, ndofs, len_t)
    for i = 0:k
        for j = 1:len_t
            u[:,j] += P[j,i+1] * coeffs[i * ndofs+1:(i+1) * ndofs]
        end
    end
    
    # Apply homogeneous Dirichlet BCs
    u = vcat(zeros(1, len_t), u, zeros(1, len_t))

    # --- Triplet containers ---
    I = Int[]
    J = Int[]
    V = ComplexF64[]

    # --- Element loop ---
    for el = 1:nelem
        element = elements[:, el]
        element_nodes = nodes[element]
        idx = mesh_to_space_map[element]

        h = abs(element_nodes[end] - element_nodes[1])

        # Quadrature points in physical space
        xg = 0.5 * h .* ξ .+ 0.5 * (element_nodes[end] + element_nodes[1])

        # Shape function matrix φ(x)
        φ = lagrange_matrix(element_nodes, xg)

        # Evaluate u(xg, t)
        u_el = u[element, :]
        ψ = φ * u_el

        # |u|² at quadrature points
        κ = zeros(ComplexF64,len_x,len_t)
        for j = 1:len_t
            κ[:,j] = abs2.(ψ[:,j]).^2
        end

        # Local integral accumulator
        int = zeros(ComplexF64, nloc, nloc)

        # Perform tensor-product quadrature in space & time
        for j = 1:len_x
            for m = 1:len_t
                for l = 1:d
                    for k = 1:d
                        int[k,l] += 0.25*h*τ*ζ[m]*ω[j]*φ[j,k]*φ[j,l]*κ[j,m]*P[m,p+1]*P[m,q+1]
                    end
                end
            end
        end

        # --- Append local contributions as triplets ---
        for jloc = 1:nloc
            gj = idx[jloc]
            if gj == 0; continue; end
            for kloc = 1:nloc
                gk = idx[kloc]
                if gk == 0; continue; end
                push!(I, gj)
                push!(J, gk)
                push!(V, int[jloc, kloc])
            end
        end
    end

    return I, J, V
end

"""
    derivative_Γ(q, k, τ, coeffs, nodes, elements, mesh_to_space_map, ndofs, side)

Compute the nonlinear derivative contribution Γ′₍q₎ for the GPE variational scheme,
returning its sparse triplet representation.

This term corresponds to the derivative of Γ with respect to time boundary
conditions (left or right) and depends on the side argument.

# Arguments
- `q::Int` : derivative order.
- `k::Int` : polynomial order in time.
- `τ::Float64` : time step.
- `coeffs::Vector{ComplexF64}` : stacked Legendre coefficients `[U₀; U₁; …; Uₖ]`.
- `nodes::Vector{Float64}` : node coordinates.
- `elements::Matrix{Int}` : FEM connectivity.
- `mesh_to_space_map::Vector{Int}` : mapping local→global.
- `ndofs::Int` : number of spatial DOFs.
- `side::Symbol` : either `:left` or `:right`.

# Returns
- `(I, J, V)::Tuple{Vector{Int}, Vector{Int}, Vector{ComplexF64}}`
"""
function derivative_Γ(q::Int, k::Int, τ::Float64,
                      coeffs::Vector{ComplexF64},
                      nodes::Vector{Float64},
                      elements::Matrix{Int},
                      mesh_to_space_map::Vector{Int},
                      ndofs::Int,
                      side::Symbol)

    nloc, nelem = size(elements)
    ξ = Quad_5.x
    ω = Quad_5.w
    len_x = length(ξ)

    # --- Precompute nonlinear field combination u(x) ---
    u = zeros(ComplexF64, ndofs)

    for s = 0:k
        for p = 0:k
            prod = imag.(coeffs[p*ndofs+1:(p+1)*ndofs]).*imag.(coeffs[s*ndofs+1:(s+1)*ndofs])+real.(coeffs[p*ndofs+1:(p+1)*ndofs]).*real.(coeffs[s*ndofs+1:(s+1)*ndofs])
            for l = 0:q
                u += binomial(q,l) * legendre_derivative_endpoint(q-l,s,τ,side=:side)* legendre_derivative_endpoint(l,p,τ,side=:side) * prod
            end
        end
    end

    # Enforce Dirichlet BCs
    u = vcat(0im, u, 0im)

    # --- Triplet containers ---
    I = Int[]
    J = Int[]
    V = ComplexF64[]

    # --- Quadrature loop ---
    for el = 1:nelem
        element = elements[:, el]
        elem_nodes = nodes[element]
        idx = mesh_to_space_map[element]

        h = abs(elem_nodes[end] - elem_nodes[1])
        xg = 0.5 * h .* ξ .+ 0.5 * (elem_nodes[end] + elem_nodes[1])

        # Lagrange basis
        φ = lagrange_matrix(elem_nodes, xg)

        # Local field projection
        u_el = u[element]
        ψ = φ * u_el

        # Local integral accumulator
        int = zeros(ComplexF64, nloc, nloc)

        for j = 1:len_x
            for l = 1:d
                for k = 1:d
                    int[k,l] += 0.5*h*ω[j]*φ[j,k]*φ[j,l]*ψ[j]
                end
            end
        end

        # Append as triplets
        for jl = 1:nloc
            gj = idx[jl]
            if gj == 0; continue; end
            for kl = 1:nloc
                gk = idx[kl]
                if gk == 0; continue; end
                push!(I, gj)
                push!(J, gk)
                push!(V, int[jl, kl])
            end
        end
    end

    return I, J, V
end