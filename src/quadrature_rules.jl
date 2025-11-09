"""
    QuadratureRule

A simple structure holding the nodes (`x`) and weights (`w`) of a quadrature rule.
"""
struct QuadratureRule
    x::Vector{Float64}
    w::Vector{Float64}
end

"""
    Quad_5 :: QuadratureRule

Gaussian quadrature rule with 5 integration points on the interval `[-1, 1]`.
"""
const Quad_5 = QuadratureRule(
        [-1/3*sqrt(5+2*sqrt(10/7)),
        -1/3*sqrt(5-2*sqrt(10/7)),
        0,
        1/3*sqrt(5-2*sqrt(10/7)),
        1/3*sqrt(5+2*sqrt(10/7))],

        [(322-13*sqrt(70))/900,
        (322+13*sqrt(70))/900,
        128/225,
        (322+13*sqrt(70))/900,
        (322-13*sqrt(70))/900])

"""
    Quad_9 :: QuadratureRule

Gaussian quadrature rule with 9 integration points on the interval `[-1, 1]`.
"""
const Quad_9 = QuadratureRule(
     [-0.9681602395076261,
        -0.8360311073266358,
        -0.6133714327005904,
        -0.3242534234038089,
        0,
        0.3242534234038089,
        0.6133714327005904,
        0.8360311073266358,
        0.9681602395076261],

     [0.0812743883615744,
        0.1806481606948574,
        0.2606106964029354,
        0.3123470770400029,
        0.3302393550012598,
        0.3123470770400029,
        0.2606106964029354,
        0.1806481606948574,
        0.0812743883615744])

"""
    legendre_gauss_lobatto(k::Integer) -> x::Vector{Float64}, w::Vector{Float64}

Compute the Legendre–Gauss–Lobatto (LGL) quadrature nodes and weights
for a polynomial of degree `k`. The resulting rule integrates exactly all
polynomials of degree ≤ `2k−3` over the interval `[-1, 1]`, and includes
both endpoints `x₁ = -1` and `xₙ = 1`.

# Arguments
- `k::Integer`: Polynomial degree (`k+1` total nodes).

# Returns
- `x::Vector{Float64}`: LGL nodes (length `k+1`, with `x[1] = -1` and `x[end] = 1`).
- `w::Vector{Float64}`: LGL quadrature weights (same length).

# Notes
The nodes are the zeros of `(1-x^2)P'_k(x)`, where `P_k` is the Legendre polynomial of degree `k`. They are found by Newton-Raphson iteration using Chebyshev-Gauss-Lobatto nodes as an initial guess.
"""
function legendre_gauss_lobatto(k::Int64)
    N = k + 1
    x = cos.(π * (0:k) ./ k) # Chebyshev–Lobatto nodes
    x_old = fill(2.0, N)
    P = zeros(N, N) # Legendre Vandermonde matrix
    
    while maximum(abs.(x .- x_old)) > 1e-12
        x_old .= x

        P[:, 1] .= 1.0
        P[:, 2]  .= x
        for j in 2:k
            P[:, j+1] = ((2j - 1) * x .* P[:, j] - (j - 1) * P[:, j-1]) / j
        end

        x .= x_old .- (x .* P[:, N] - P[:, k]) ./ (N * P[:, N])
    end

    w = 2.0 ./ (k * N * P[:, N].^2)

    return x, w
end

"""
    legendre_gauss_radau(k::Int64)

Compute the Legendre–Gauss–Radau (LGR) quadrature nodes and weights for
a polynomial of degree `k`. The resulting rule integrates exactly all
polynomials of degree ≤ `2k−1` over the interval `[-1, 1]`, with one
endpoint fixed at `x₁ = -1`.

# Arguments
- `k::Integer`: polynomial degree (number of free nodes = `k`).

# Returns
- `x::Vector{Float64}`: LGR nodes (length `k+1`, first node = -1).
- `w::Vector{Float64}`: LGR quadrature weights (same length).

# Notes
The nodes are obtained via Newton–Raphson iteration on the recurrence
relation for Legendre polynomials. Initial guesses are Chebyshev–Gauss–Radau
nodes, which ensure rapid convergence.
"""
function legendre_gauss_radau(k::Int64)
    N = k + 1
    x = -cos.(2π * (0:k) ./ (2k + 1)) # Chebyshev–Radau nodes
    x_old = fill(2.0, N) 

    P = zeros(N, N + 1)
    free = 2:N                  

    while maximum(abs.(x .- x_old)) > 1e-12
        x_old .= x

        P[1, :] .= (-1).^(0:N)
        P[free, 1] .= 1.0
        P[free, 2]  .= x[free]
        for j in 2:N
            P[free, j+1] = ((2j - 1) * x[free] .* P[free, j] - (j - 1) * P[free, j-1]) / j
        end

        x[free] = x_old[free] - ((1 .- x_old[free]) / N) .* (P[free, N] + P[free, N+1]) ./ (P[free, N] - P[free, N+1])
        x[1] = -1.0                
    end

    w = zeros(N)
    w[1] = 2 / N^2
    w[free] = (1 .- x[free]) ./ (N * P[free, N]).^2

    return x, w
end
