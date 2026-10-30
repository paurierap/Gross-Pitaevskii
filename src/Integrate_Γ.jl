"""
    NonlinearCache(nodes, elements, map, ndofs, k, r, τ)

Precompute space/time quadrature, basis values and sparse triplet indices.
All time evaluations use [0,τ], avoiding cancellation at large physical times.
"""
struct NonlinearCache
    basis::Vector{Matrix{Float64}}
    weights::Vector{Vector{Float64}}
    dofs::Matrix{Int}
    P::Matrix{Float64}
    wt::Vector{Float64}
    left::Matrix{Float64}
    right::Matrix{Float64}
    rows::Vector{Tuple{Int,Int,Symbol}}
    I::Vector{Int}
    J::Vector{Int}
    ndofs::Int
    k::Int
end

function NonlinearCache(nodes, elements, map, ndofs, k, r, τ)
    nloc, nelem = size(elements)
    ξ, ω = gauss_legendre(2*(nloc-1)+1) # Exact for degree 4d.
    η, ζ = gauss_legendre(2k+1)        # Exact for degree 4k.
    P = legendre_matrix(k, [0.0, τ], τ/2 .* (η .+ 1))
    basis = Matrix{Float64}[]
    weights = Vector{Float64}[]
    for el in 1:nelem
        xn = nodes[elements[:,el]]
        h = xn[end]-xn[1]
        # Reference coordinates keep evaluation well-conditioned on fine meshes.
        ref = 2 .* (xn .- xn[1]) ./ h .- 1
        push!(basis, lagrange_matrix(ref, ξ))
        push!(weights, h/2 .* ω)
    end
    D0, D1 = (r-1)÷2, r÷2
    left = [legendre_derivative_endpoint(q, s, τ; side=:left) for s in 0:k, q in 0:max(D1-1,0)]
    right = [legendre_derivative_endpoint(q, s, τ; side=:right) for s in 0:k, q in 0:max(D1-1,0)]
    rows = [(p, -1, :interior) for p in 0:(k-r)]
    append!(rows, [(k-r+1+q, q-1, :left) for q in 1:D0])
    append!(rows, [(k-r+1+D0+q, q-1, :right) for q in 1:D1])
    dofs = map[elements]
    I, J = Int[], Int[]
    for el in 1:nelem, (row, _, _) in rows, col in 0:k, b in 1:nloc, a in 1:nloc
        ia, ib = dofs[a,el], dofs[b,el]
        (ia == 0 || ib == 0) && continue
        push!(I, row*ndofs+ia)
        push!(J, col*ndofs+ib)
    end
    return NonlinearCache(basis, weights, dofs, P, τ/2 .* ζ, left, right, rows, I, J, ndofs, k)
end

# Reconstruct the field before taking products: interpolation of |nodal u|²
# would omit cross terms between different spatial basis functions.
function element_coefficients(coeffs, cache, el)
    C = reshape(coeffs, cache.ndofs, cache.k+1)
    local_coeffs = zeros(ComplexF64, size(cache.dofs,1), cache.k+1)
    for a in axes(local_coeffs,1)
        i = cache.dofs[a,el]
        i == 0 && continue
        local_coeffs[a,:] .= @view C[i,:]
    end
    return cache.basis[el] * local_coeffs
end

function density_derivatives(U, endpoints)
    derivatives = U * endpoints
    ρ = zeros(size(derivatives))
    for q in 0:(size(ρ,2)-1), α in 0:q
        @views ρ[:,q+1] .+= binomial(q,α) .* real.(conj.(derivatives[:,α+1]) .* derivatives[:,q-α+1])
    end
    return ρ
end

"""
    nonlinear_matrix(coeffs, cache)

Assemble the nonlinear block matrix Γ(u) without β.
Each element's space-time field is reconstructed only once per iteration.
"""
function nonlinear_matrix(coeffs, cache::NonlinearCache)
    V = Vector{Float64}(undef, length(cache.I))
    index = 0
    nloc, nelem = size(cache.dofs)
    has_endpoints = any(row[2] ≥ 0 for row in cache.rows)
    for el in 1:nelem
        φ, wx = cache.basis[el], cache.weights[el]
        U = element_coefficients(coeffs, cache, el)
        density = abs2.(U * transpose(cache.P))
        ρleft = has_endpoints ? density_derivatives(U, cache.left) : nothing
        ρright = has_endpoints ? density_derivatives(U, cache.right) : nothing
        for (row, q, side) in cache.rows, col in 0:cache.k
            if side == :interior
                moment = density * (cache.wt .* cache.P[:,row+1] .* cache.P[:,col+1])
            else
                ρ = side == :left ? ρleft : ρright
                E = side == :left ? cache.left : cache.right
                moment = zeros(length(wx))
                for α in 0:q
                    @views moment .+= binomial(q,α) * E[col+1,α+1] .* ρ[:,q-α+1]
                end
            end
            local_matrix = transpose(φ) * ((wx .* moment) .* φ)
            for b in 1:nloc, a in 1:nloc
                (cache.dofs[a,el] == 0 || cache.dofs[b,el] == 0) && continue
                index += 1
                V[index] = local_matrix[a,b]
            end
        end
    end
    n = (cache.k+1)*cache.ndofs
    return sparse(cache.I, cache.J, V, n, n)
end

"""
    integrate_Γ(coeffs, p, q, k, time_interval, nodes, elements, map, ndofs)

Return spatial triplets for a single variational block Γᵖᑫ, without β.
The time-stepper uses `nonlinear_matrix` to assemble all blocks together.
"""
function integrate_Γ(coeffs, p, q, k, time_interval, nodes, elements, map, ndofs)
    τ = time_interval[2]-time_interval[1]
    cache = NonlinearCache(nodes, elements, map, ndofs, k, 1, τ)
    I, J, V = Int[], Int[], Float64[]
    for el in axes(elements,2)
        φ = cache.basis[el]
        U = element_coefficients(coeffs, cache, el)
        density = abs2.(U * transpose(cache.P))
        moment = density * (cache.wt .* cache.P[:,p+1] .* cache.P[:,q+1])
        local_matrix = transpose(φ) * ((cache.weights[el] .* moment) .* φ)
        append_spatial_triplets!(I, J, V, local_matrix, cache.dofs[:,el])
    end
    return I, J, V
end

"""
    derivative_Γ(q, k, τ, coeffs, nodes, elements, map, ndofs; side=:right)

Return spatial triplets for the q-th time derivative of the density matrix.
Products are evaluated at quadrature points using Leibniz's rule.
"""
function derivative_Γ(q, k, τ, coeffs, nodes, elements, map, ndofs; side=:right)
    @assert q ≥ 0
    cache = NonlinearCache(nodes, elements, map, ndofs, k, 1, τ)
    E = [legendre_derivative_endpoint(α, s, τ; side) for s in 0:k, α in 0:q]
    I, J, V = Int[], Int[], Float64[]
    for el in axes(elements,2)
        φ = cache.basis[el]
        ρ = density_derivatives(element_coefficients(coeffs, cache, el), E)
        local_matrix = transpose(φ) * ((cache.weights[el] .* ρ[:,q+1]) .* φ)
        append_spatial_triplets!(I, J, V, local_matrix, cache.dofs[:,el])
    end
    return I, J, V
end

function append_spatial_triplets!(I, J, V, local_matrix, dofs)
    for b in eachindex(dofs), a in eachindex(dofs)
        (dofs[a] == 0 || dofs[b] == 0) && continue
        push!(I, dofs[a]); push!(J, dofs[b]); push!(V, local_matrix[a,b])
    end
end
