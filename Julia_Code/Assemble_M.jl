function Assemble_M(Nodes::Array{Float64,1},Simplices::Array{Int64,2},Mesh2Space::Array{Int64,1},SpaceSize::Int64)

    M = spzeros(SpaceSize,SpaceSize)

    d = size(Simplices,1)

    N = size(Simplices,2)

    for i = 1:N

        Simplex = Simplices[:,i]

        Element_Nodes = Nodes[Simplex]

        loc = Assemble_Element_M(Element_Nodes) # Calculate Element matrix, integrated exactly with Gaussian quadrature.

        idx = Mesh2Space[Simplex]

        for j = 1:d

            for k = 1:d

                if (idx[j]!=0) & (idx[k]!=0)

                    M[idx[j],idx[k]] =  M[idx[j],idx[k]]+loc[j,k]

                end

            end

        end

    end

    return M

end
