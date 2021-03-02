function energy_integration(β,u,Nodes,Simplices)

    # Assume Dirichlet BC's.

    u = [0;u;0]

    energy = 0

    d = size(Simplices,1)

    N = size(Simplices,2)

    Id = Matrix(I,d,d) # Identity matrix of size dxd.

    # 5-point Gaussian:

    ξ = Quad_5["X"]
    ω = Quad_5["ω"]

    for i = 1:N

        Simplex = Simplices[:,i]

        Element_Nodes = Nodes[Simplex]
        u_τ = u[Simplex]

        h = abs(Element_Nodes[end]-Element_Nodes[1])
        x = h*0.5*ξ.+0.5*(Element_Nodes[end]+Element_Nodes[1])

        int = 0

        V = ones(d,d)
        X = ones(length(x),d)

        for j = 2:d

            V[:,j] = V[:,j-1].*Element_Nodes
            X[:,j] = X[:,j-1].*x

        end

        φ = X*(V\u_τ) # A little bit of algebra required here.

        L = LagDerivative(Element_Nodes,x)

        δφ = L*u_τ

        for j = 1:length(x)

            int += ω[j]*(abs(δφ[j])^2+0.5*β*abs(φ[j])^4)

        end

        energy += 0.5*h*int

    end

    return energy

end
