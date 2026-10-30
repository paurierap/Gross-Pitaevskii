using Test, LinearAlgebra, SparseArrays
include("../src/LoadCode.jl")

"Analytic two-soliton solution of iuₜ = -uₓₓ - 2|u|²u."
function exact_solution(x, t)
    den = -128 .* cos.(12 .* t) .+ 4 .* exp.(-6 .* x) .+ 16 .* exp.(6 .* x) .+
          81 .* exp.(-2 .* x) .+ 64 .* exp.(2 .* x)
    num = 8 .* exp.(4im .* t) .* (9 .* exp.(-4 .* x) .+ 16 .* exp.(4 .* x)) .-
          32 .* exp.(16im .* t) .* (4 .* exp.(-2 .* x) .+ 9 .* exp.(2 .* x))
    return num ./ den
end

function l2_error(nodes, elements, u, t)
    u = vcat(0im, u, 0im)
    ξ, ω = gauss_legendre(12)
    error = 0.0
    for el in axes(elements, 2)
        element_nodes = nodes[elements[:, el]]
        h = element_nodes[end]-element_nodes[1]
        xq = 0.5 .* h .* ξ .+ 0.5 .* (element_nodes[end]+element_nodes[1])
        φ = lagrange_matrix(2 .* (element_nodes .- element_nodes[1]) ./ h .- 1, ξ)
        uq = φ * u[elements[:, el]]
        error += 0.5 * h * sum(ω .* abs2.(uq .- exact_solution(xq, t)))
    end
    return sqrt(error)
end

@testset "Legendre basis is orthonormal" begin
    τ = 0.3
    ξ, ω = gauss_legendre(8)
    t = 0.5 .* τ .* (ξ .+ 1)
    P = legendre_matrix(4, [0.0, τ], t)
    @test P' * Diagonal(0.5 .* τ .* ω) * P ≈ I(5) atol = 1e-13
end

@testset "Nonlinear density is quadratic" begin
    # A stationary P1 hat on [0,2]: ∫ |u|² φ² dx = ∫ φ⁴ dx = 2/5.
    nodes = [0.0, 1.0, 2.0]
    elements = [1 2; 2 3]
    map = [0, 1, 0]
    coeffs = ComplexF64[1, 0]
    I, J, V = integrate_Γ(coeffs, 0, 0, 1, [0.0, 1.0], nodes, elements, map, 1)
    @test sparse(I, J, V, 1, 1)[1, 1] ≈ 2/5
end

@testset "Full nonlinear matrix matches individual blocks" begin
    nodes = collect(0.0:0.5:2.0)
    elements = [1 2 3 4; 2 3 4 5]
    map = [0; 1; 2; 3; 0]
    ndofs, k, τ = 3, 3, 0.2
    coeffs = ComplexF64[sin(i) + im*cos(2i) for i in 1:(k+1)*ndofs]
    cache = NonlinearCache(nodes, elements, map, ndofs, k, 1, τ)
    Γ = nonlinear_matrix(coeffs, cache)
    for p in 0:(k-1), q in 0:k
        I, J, V = integrate_Γ(coeffs, p, q, k, [0.0, τ], nodes, elements, map, ndofs)
        @test Γ[block_inds(p, ndofs), block_inds(q, ndofs)] ≈
              sparse(I, J, V, ndofs, ndofs) atol = 1e-13
    end
end

@testset "Endpoint collocation blocks use the correct side" begin
    nodes = collect(0.0:0.5:2.0)
    elements = [1 2 3 4; 2 3 4 5]
    map = [0; 1; 2; 3; 0]
    ndofs, k, τ = 3, 4, 0.2
    coeffs = ComplexF64[sin(i) + im*cos(2i) for i in 1:(k+1)*ndofs]
    for r in (2, 4)
        cache = NonlinearCache(nodes, elements, map, ndofs, k, r, τ)
        Γ = nonlinear_matrix(coeffs, cache)
        D0, D1 = (r-1)÷2, r÷2
        for (rows, side) in ((1:D0, :left), (1:D1, :right))
            offset = side == :left ? 0 : D0
            for q in rows
                row = k-r+1+offset+q
                for col in 0:k
                    expected = spzeros(ndofs, ndofs)
                    for α in 0:(q-1)
                        I, J, V = derivative_Γ(q-1-α, k, τ, coeffs,
                                               nodes, elements, map, ndofs; side)
                        expected += binomial(q-1, α) *
                                    legendre_derivative_endpoint(α, col, τ; side) *
                                    sparse(I, J, V, ndofs, ndofs)
                    end
                    @test Γ[block_inds(row, ndofs), block_inds(col, ndofs)] ≈ expected
                end
            end
        end
    end
end

@testset "Mesh endpoints and requested degree" begin
    for degree in 1:4, nx in (7, 8), spacing in (:uniform, :chebyshev)
        nodes, elements, t, intervals = generate_mesh_1D(-20.0, 20.0, nx, 0.0, 1.0, 4; d=degree, spacing)
        @test length(nodes) == nx * degree + 1
        @test first(nodes) == -20.0
        @test last(nodes) == 20.0
        @test all(diff(nodes) .> 0)
        @test maximum(elements) == length(nodes)
    end
end

@testset "Analytic two-soliton solution" begin
    T, d, k, r = 0.2, 2, 3, 1
    nodes, elements, _, intervals = generate_mesh_1D(-10.0, 10.0, 128, 0.0, T, 32; d=d)
    ndofs = length(nodes)-2
    map = [0; collect(1:ndofs); 0]
    S = assemble_stiffness_matrix(nodes, elements, map, ndofs)
    M = assemble_mass_matrix(nodes, elements, map, ndofs)
    u0 = exact_solution(nodes[2:end-1], 0.0)
    u, mass, energy = Gr(M, S, intervals, u0, nodes, elements, map, ndofs, k, r, -2.0; tol=1e-9)
    @test l2_error(nodes, elements, u[:, 1], T) < 0.5
    @test abs(mass[end]-mass[1]) < 1e-5
    @test abs(energy[end]-energy[1]) < 1e-4
end
