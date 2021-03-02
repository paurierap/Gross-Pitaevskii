function Derivative_M_tau(q::Int64,k::Int64,τ::Float64,coeffs,Nodes::Array{Float64,1},Simplices::Array{Int64,2},Mesh2Space::Array{Int64,1},SpaceSize::Int64,side::String)

    M_tau = spzeros(ComplexF64,SpaceSize,SpaceSize)

    d = size(Simplices,1)

    N = size(Simplices,2)

    Id = Matrix(I,d,d) # Identity matrix of size dxd.

    u = zeros(ComplexF64,SpaceSize,1)

    #ψ = [zeros(ComplexF64,1,k+1);reshape(coeffs,SpaceSize,k+1);zeros(ComplexF64,1,k+1)]

    # Create derivative pattern:

    for s = 0:k

        for p = 0:k

            prod = coeffs[p*SpaceSize+1:(p+1)*SpaceSize].*coeffs[s*SpaceSize+1:(s+1)*SpaceSize]

            for l = 0:q

                u += binomial(q,l)*LegDerivative(q-l,s,τ,side)*LegDerivative(l,p,τ,side)*prod

            end

        end

    end

    u = [0;u;0]

    # 5-point Gaussian:

    ξ = Quad_5["X"]
    ω = Quad_5["ω"]

    for i = 1:N

        Simplex = Simplices[:,i]

        Element_Nodes = Nodes[Simplex]
        u_τ = u[Simplex]

        h = abs(Element_Nodes[end]-Element_Nodes[1])
        x = h*0.5*ξ.+0.5*(Element_Nodes[end]+Element_Nodes[1])

        int = zeros(ComplexF64,d,d)

        V = ones(d,d)
        X = ones(length(x),d)

        for j = 2:d

            V[:,j] = V[:,j-1].*Element_Nodes
            X[:,j] = X[:,j-1].*x

        end

        #=cross = 0

        for s = 0:k

            for p = 0:k

                prod = ψ[Simplex[1],s+1]*ψ[Simplex[2],p+1]+ψ[Simplex[2],s+1]*ψ[Simplex[1],p+1]

                for l = 0:q

                    cross += binomial(q,l)*LegDerivative(q-l,s,τ,side)*LegDerivative(l,p,τ,side)*prod

                end

            end

        end

        u_τ = u[Simplex[1]]*(1/h*(b.-x)).^2 + cross/h.*(b.-x)/h.*(x.-a) + u[Simplex[2]]*(1/h*(x.-a)).^2

        =#

        φ = X*(V\Id) # Basis functions.

        ψ = X*(V\u_τ) # Approximated function.

        for j = 1:length(x)

            for l = 1:d

                for k = 1:d

                    int[k,l] += ω[j]*φ[j,k]*φ[j,l]*ψ[j]

                end

            end

        end

        int = 0.5*h*int

        idx = Mesh2Space[Simplex];

        for j = 1:d

            for k = 1:d

                if ((idx[j]!=0) & (idx[k]!=0))

                    M_tau[idx[j],idx[k]] =  M_tau[idx[j],idx[k]]+int[j,k]

                end

            end

        end

    end

    return M_tau

end
