"""
    legendre_matrix(k::Int64,In::Vector{Float64},t::Vector{Float64})

Evaluate the first k+1 Legendre polynomials (including `L_0`) for a given interval at a set of evaluation points.

Given an interval `In = [t_{n-1}, t_n]`, this function evaluates each Legendre Polynomial `L_i(x)` for `i=0,...,k` at the positions `t`. The resulting matrix `L` has size `(length(t), k+1)` such that `L[i, j] = L_{j-1}(t[i])`.

# Arguments
- `k::Int64`: index of the last Legendre polynomial
- `In::Vector{Float64}`: time interval considered.
- `t::Vector{Float64}`: time points at which the polynomials are evaluated 

# Returns
L::{Matrix, Float64}, where each column represents the evaluation of the corresponding Legendre Polynomial
"""
function legendre_matrix(k::Int64, time_interval::Vector{Float64}, t::Vector{Float64})
    @assert k ≥ 0 "Polynomial order k must be non-negative."
    L = zeros(length(t), k+1);
    τ = abs(time_interval[2] - time_interval[1]);

    L[:,1] .= 1/sqrt(τ); 
    if k > 0
        L[:,2] = sqrt(3/τ)*(2*(t .- time_interval[1])/τ .-1);

        for j = 3:(k+1)
            L[:,j] = sqrt((2*j-1)/τ)/(j-1) * ((2*j-3)*sqrt(τ/(2*j-3)) * (2*(t .- time_interval[1])/τ .-1) .* L[:,j-1] .- sqrt(τ/(2*j-5))*(j-2) * L[:,j-2]);
        end
    end

    return L
end

"""
    lagrange_matrix(nodes::Vector{Float64}, x_eval::Vector{Float64})

Given interpolation nodes `nodes = [x₁, x₂, …, xₙ]`, evaluates each Lagrange basis polynomial `L_j(x)` at all positions
in `x_eval`. The resulting matrix `L` has size `(length(x_eval), length(nodes))`
such that `L[k, j] = L_j(x_eval[k])`.

# Arguments
- `nodes::Vector{Float64}`: Interpolation nodes.
- `x_eval::Vector{Float64}`: Points where the polynomials are to be evaluated.

# Returns
- `L::Matrix{Float64}`: Lagrange polynomial matrix of size `(length(x_eval), length(nodes))`.
"""
function lagrange_matrix(nodes::Vector{Float64}, x_eval::Vector{Float64})
    nloc = length(nodes)

    # Build Vandermonde matrix for polynomial basis
    V = ones(nloc, nloc)
    X = ones(length(x_eval), nloc)
    for j = 2:nloc
        V[:, j] = nodes .^ (j-1)
        X[:, j] = x_eval .^ (j-1)
    end

    # Compute Lagrange basis functions at quadrature points: φ(q,:) = X * inv(V)
    return X * (V \ I(nloc))
end

"""
    lagrange_derivative_matrix(nodes::Vector{Float64}, x_eval::Vector{Float64})

Compute the derivatives of the Lagrange basis polynomials at a set of evaluation points.

Given interpolation nodes `nodes = [x₁, x₂, …, xₙ]`, this function computes
the derivative of each Lagrange basis polynomial `L_j(x)` at all positions
in `x_eval`. The resulting matrix `D` has size `(length(x_eval), length(nodes))`
such that `D[k, j] = L_j'(x_eval[k])`.

# Arguments
- `nodes::Vector{Float64}`: Interpolation nodes.
- `x_eval::Vector{Float64}`: Points where the derivatives are to be evaluated.

# Returns
- `D::Matrix{Float64}`: Derivative matrix of size `(length(x_eval), length(nodes))`.
"""
function lagrange_derivative_matrix(nodes::Vector{Float64}, x_eval::Vector{Float64})
    n_nodes = length(nodes)
    n_points = length(x_eval)
    D = zeros(n_points, n_nodes)
    temp = ones(n_points)

    for j = 1:n_nodes
        for i = 1:n_nodes
            fill!(temp, 1.0)
            for m in 1:n_nodes
                if m != i && m != j
                    temp .*= (x_eval .- nodes[m]) / (nodes[j] - nodes[m])
                end
            end
            if i != j
                D[:, j] .+= (1 / (nodes[j] - nodes[i])) .* temp
            end
        end
    end

    return D
end

"""
    legendre_derivative_endpoint(q::Int64, s::Int64, τ::Float64; side::Symbol = :right)

Compute the value of the q-th derivative of the s-th Legendre polynomial,
scaled to an interval of length `τ`, evaluated at either the left or right endpoint.

This returns:
    P_s^(q)(x) evaluated at either ±τ/2

where `side = :left` or `:right` selects which endpoint to use.
The scaling from `[-1, 1]` to `[-τ/2, τ/2]` is included through the factor `(2/τ)^q`.

# Arguments
- `q::Integer`: Derivative order.
- `s::Integer`: Polynomial degree.
- `τ::Real`: Interval length.
- `side::Symbol`: Either `:left` or `:right` (default = `:right`).

# Returns
- `Float64`: Value of the scaled derivative at the chosen endpoint.
"""
function legendre_derivative_endpoint(q::Int64, s::Int64, τ::Float64; side::Symbol = :right)
    @assert s ≥ 0 "Polynomial order s must be non-negative."
    @assert q ≥ 0 "Derivative order q must be non-negative."
    @assert s ≥ q "Derivative order q cannot exceed polynomial degree s."
    @assert τ > 0 "Interval length τ must be positive."
    @assert side in (:left, :right) "Side must be :left or :right."

    prefactor = sqrt((2s + 1) / τ) * (2 / τ)^q * derivative_factor(q) * binomial(s + q, s - q)

    return side === :left ? (-1)^(s - q) * prefactor : prefactor
end

"""
    derivative_factor(q::Int64)
Compute the product of odd integers up to `2q - 1`:
∏_{k=1}^{q} (2k - 1)
Used in analytic expressions for Legendre polynomial derivatives.
"""
@inline derivative_factor(q::Int64) = prod(1:2:(2q - 1))

"""
    is_coupled(i::Int64, j::Int64)

Return `true` if the inner product ∫ Pᵢ'(x) Pⱼ(x) dx is nonzero.

This happens when the derivative of the i-th Legendre polynomial couples
with the j-th polynomial, i.e. when `i > j` and `i - j` is odd.
"""
@inline function is_coupled(i::Int64, j::Int64)
    return i > j && isodd(i - j)
end