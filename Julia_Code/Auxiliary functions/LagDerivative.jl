function LagDerivative(Nodes::Array{Float64,1},x::Array{Float64,1})

    # Calculates Lagrange Polynomials 1st derivative for Pd-FEM.

    d = length(Nodes)
    q = length(x)

    L = zeros(q,d)

    for j = 1:d

        for i = 1:d

            l = ones(q,1)

            for m = 1:d

                if (m != i) & (m != j)

                    l .*= (x .- Nodes[m])/(Nodes[j] - Nodes[m])

                end

            end

            if i != j

                L[:,j] += 1/(Nodes[j]-Nodes[i])*l

            end

        end

    end

  return L

end
