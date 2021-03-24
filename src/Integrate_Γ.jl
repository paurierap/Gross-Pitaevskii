function Integrate_Γ(coeffs::Array{ComplexF64,1},p::Int64,q::Int64,k::Int64,In::Array{Float64,1},Nodes::Array{Float64,1},Simplices::Array{Int64,2},Mesh2Space::Array{Int64,1},SpaceSize::Int64)

    Γ = spzeros(SpaceSize,SpaceSize)

    d,N = size(Simplices)

    Id = Matrix(I,d,d) # Identity matrix of size dxd.

    int = zeros(ComplexF64,d,d)

    V = zeros(ComplexF64,d,d)

    # 5-point Gaussian space quadrature:

    ξ = Quad_5["X"]
    ω = Quad_5["ω"]

    len_x = length(ξ)

    X = zeros(ComplexF64,len_x,d)

    # k+1-point Gauss-Legendre time quadrature:

    η,ζ = GLnw(k+1)

    t_n1 = In[1]; t_n = In[2]
    τ = abs(t_n-t_n1)

    t = τ*0.5*η.+0.5*(t_n1+t_n)
    len_t = length(t)

    u = zeros(ComplexF64,SpaceSize,len_t)

    κ = zeros(ComplexF64,len_x,len_t)

    P = LegPols(k,In,t)

    for i = 0:k
        for j = 1:len_t
            u[:,j] += P[j,i+1]*coeffs[i*SpaceSize+1:(i+1)*SpaceSize]
        end
    end

    u = [zeros(1,len_t);u;zeros(1,len_t)]

    for i = 1:N

        Simplex = Simplices[:,i]

        Element_Nodes = Nodes[Simplex]
        u_τ = u[Simplex,:]

        h = abs(Element_Nodes[end]-Element_Nodes[1])
        x = h*0.5*ξ.+0.5*(Element_Nodes[end]+Element_Nodes[1])

        fill!(int,0)
        fill!(V,1)
        fill!(X,1)
        fill!(κ,0)

        for j = 2:d
            V[:,j] = V[:,j-1].*Element_Nodes
            X[:,j] = X[:,j-1].*x
        end

        φ = X*(V\Id) # Basis functions.

        ψ = X*(V\u_τ) # Approximated function.

        for j = 1:len_t
            κ[:,j] = abs.(ψ[:,j]).^2
        end

        for j = 1:len_x
            for m = 1:len_t
                for l = 1:d
                    for k = 1:d
                        int[k,l] += 0.25*h*τ*ζ[m]*ω[j]*φ[j,k]*φ[j,l]*κ[j,m]*P[m,p+1]*P[m,q+1]
                    end
                end
            end
        end

        idx = Mesh2Space[Simplex];

        for j = 1:d
            for k = 1:d
                if idx[j] != 0 && idx[k] != 0
                    Γ[idx[j],idx[k]] += int[j,k]
                end
            end
        end

    end

    return Γ
end
