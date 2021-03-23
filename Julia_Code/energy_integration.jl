function energy_integration(β::Float64,u::Array{ComplexF64,1},Nodes::Array{Float64,1},Simplices::Array{Int64,2})

    # Assume Dirichlet BC's.

    u = [0;u;0]

    energy = 0.

    d,N = size(Simplices)

    Id = Matrix(I,d,d) # Identity matrix of size dxd.

    V = zeros(d,d)

    # 5-point Gaussian:

    ξ = Quad_5["X"]
    ω = Quad_5["ω"]

    len_x = length(ξ)

    X = zeros(len_x,d)

    for i = 1:N

        Simplex = Simplices[:,i]

        Element_Nodes = Nodes[Simplex]
        u_τ = u[Simplex]

        h = abs(Element_Nodes[end]-Element_Nodes[1])
        x = h*0.5*ξ.+0.5*(Element_Nodes[end]+Element_Nodes[1])

        fill!(V,1)
        fill!(X,1)

        for j = 2:d
            V[:,j] = V[:,j-1].*Element_Nodes
            X[:,j] = X[:,j-1].*x
        end

        φ = X*(V\u_τ) # A little bit of algebra required here.

        L = LagDerivative(Element_Nodes,x);

        δφ = L*u_τ

        int = 0.

        for j = 1:len_x
            int += 0.5*h*ω[j]*(abs(δφ[j])^2+0.5*β*abs(φ[j])^4)
        end

        energy += int

    end

    return energy

end
