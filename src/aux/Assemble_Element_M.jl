function Assemble_Element_M(Nodes::Array{Float64,1})

    d = length(Nodes)

    M = zeros(d,d)

    ξ,ω = GLnw(d) # Gauss-Lobatto quadrature. It provides d+1 points and integrates exactly polynomials of order 2d-1.

    h = abs(Nodes[end]-Nodes[1])

    x = h*0.5*ξ.+0.5*(Nodes[end]+Nodes[1])

    V = ones(d,d)
    X = ones(length(x),d)
    Id = Matrix(I,d,d) # Identity matrix of size dxd.

    for j = 2:d
        V[:,j] = V[:,j-1].*Nodes
        X[:,j] = X[:,j-1].*x
    end

    φ = X*(V\Id) # A little bit of algebra required here.

    for i = 1:d
        for j = i:d
            for k = 1:(d+1)
                M[i,j] += 0.5*h*ω[k]*φ[k,i]*φ[k,j]
            end
            M[j,i] = M[i,j]
        end
    end

    return M
end
