function Assemble_M_tau(u0,Nodes::Array{Float64,1},Simplices::Array{Int64,2},Mesh2Space::Array{Int64,1},SpaceSize::Int64,N_of_Simplices::Int64)

    M_tau = spzeros(SpaceSize,SpaceSize)

    u = [0;u0;0]

    # 5-point Gaussian:

    ξ = Quad_5["X"]
    ω = Quad_5["ω"]

    for i = 1:N_of_Simplices

        Simplex = Simplices[:,i]

        Vertexes = Nodes[Simplex]

        a = Vertexes[1]
        b = Vertexes[2]

        h = abs(b-a)

        int = zeros(2,2)

        V = inv([ones(2,1) Vertexes])

        x = h*0.5*ξ.+0.5*(a+b)

        u_tau = abs.(u[Simplex[1]]/h*(b.-x)+u[Simplex[2]]/h*(x.-a)).^2

        φ = V'*[ones(1,length(x));x']

        for j = 1:length(x)
            for l = 1:2
                for k = 1:2
                    int[k,l] += ω[j]*φ[k,j]*φ[l,j]*u_tau[j]
                end
            end
        end

        int = 0.5*h*int

        idx = Mesh2Space[Simplex];

        for j = 1:2

            for i = 1:2

                if ((idx[i]!=0) & (idx[j]!=0))

                    M_tau[idx[i],idx[j]] =  M_tau[idx[i],idx[j]]+int[i,j]

                end

            end

        end

    end

    return M_tau

end
