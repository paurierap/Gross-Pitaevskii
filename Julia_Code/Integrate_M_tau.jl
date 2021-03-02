function Integrate_M_tau(coeffs,p,q,k,In,Nodes::Array{Float64,1},Simplices::Array{Int64,2},Mesh2Space::Array{Int64,1},SpaceSize::Int64)

    M_tau = spzeros(SpaceSize,SpaceSize)

    d = size(Simplices,1)

    N = size(Simplices,2)

    Id = Matrix(I,d,d) # Identity matrix of size dxd.

    # 5-point Gaussian spatial quadrature:

    ξ = Quad_5["X"]
    ω = Quad_5["ω"]

    # k+1-point Gauss-Legendre time quadrature:

    η,ζ = GLnw(k+1);

    t_n1 = In[1]
    t_n = In[2]
    τ = abs(t_n-t_n1)

    t = τ*0.5*η.+0.5*(t_n1+t_n)

    u = zeros(ComplexF64,SpaceSize,length(t))

    P = LegPols(k,In,t);

    for i = 0:k

        for j = 1:length(t)

            u[:,j] += P[j,i+1]*coeffs[i*SpaceSize+1:(i+1)*SpaceSize]

        end

    end

    u = [zeros(1,length(t));u;zeros(1,length(t))]

    for i = 1:N

        Simplex = Simplices[:,i]

        Element_Nodes = Nodes[Simplex]
        u_τ = u[Simplex,:]

        h = abs(Element_Nodes[end]-Element_Nodes[1])
        x = h*0.5*ξ.+0.5*(Element_Nodes[end]+Element_Nodes[1])

        int = zeros(ComplexF64,d,d)

        V = ones(d,d)
        X = ones(length(x),d)

        for j = 2:d

            V[:,j] = V[:,j-1].*Element_Nodes
            X[:,j] = X[:,j-1].*x

        end

        φ = X*(V\Id) # Basis functions.

        ψ = X*(V\u_τ) # Approximated function.

        κ = zeros(length(x),length(t))

        for j = 1:length(t)

            κ[:,j] = abs.(ψ[:,j]).^2

        end

        for j = 1:length(x)

            for m = 1:length(t)

                for l = 1:d

                    for k = 1:d

                        int[k,l] += ζ[m]*ω[j]*φ[j,k]*φ[j,l]*κ[j,m]*P[m,p+1]*P[m,q+1]

                    end

                end

            end

        end

        int = 0.5*h*0.5*τ*int;

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
