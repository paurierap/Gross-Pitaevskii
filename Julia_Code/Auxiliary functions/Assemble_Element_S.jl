function Assemble_Element_S(Nodes::Array{Float64,1})

    d = length(Nodes)

    S = zeros(d,d)

    ξ,ω = GLnw(d-1) # Gauss-Lobatto quadrature. It provides d points and integrates exactly polynomials of order 2d-3.

    h = abs(Nodes[end]-Nodes[1])

    x = h*0.5*ξ.+0.5*(Nodes[end]+Nodes[1])

    L = LagDerivative(Nodes,x) # L has size (d,d).

    for i = 1:d

        for j = i:d

            for k = 1:d

                S[i,j] += ω[k]*L[k,i]*L[k,j]

            end

            S[i,j] *= 0.5*h

            S[j,i] = S[i,j]

        end

    end

    return S

end
