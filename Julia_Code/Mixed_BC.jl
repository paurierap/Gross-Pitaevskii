function Mixed_BC(boundary_left::Float64,boundary_right::Float64,Nodes::Array{Float64,1},kind::String)

    N = length(Nodes);

    Mesh2Space = zeros(Int64,N);
    fill!(Mesh2Space,0);

    idx = 1;

    if kind == "Dirichlet"

        for j = 1:N
            if (Nodes[j] == boundary_left) | (Nodes[j] == boundary_right)

            else
                Mesh2Space[j] = idx;
                idx = idx+1
            end
        end

    elseif kind == "Neumann"

        for j = 1:N
            Mesh2Space[j] = idx;
            idx = idx+1
        end

    elseif kind == "Robin"

        for j = 1:N
            Mesh2Space[j] = idx;
            idx = idx+1;
        end

    end

    SpaceSize = idx-1

    return Mesh2Space, SpaceSize

end
