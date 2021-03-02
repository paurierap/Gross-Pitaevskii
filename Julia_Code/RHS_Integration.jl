function RHS_Integration(a,b,k,D0,D1,f,Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,kind)

    # (k+1)-point Gauss-Radau or Gauss-Lobatto quadratures.

    if kind == "dG"

        ξ,ω = GRnw(k);

    elseif kind == "cGP"

        ξ,ω = GLnw(k+1);

    elseif kind == "G"

        ξ = zglj(k+1,D0+1,D1+1);
        ω = wglj(chi,D0+1,D1+1);

    end


    t = (b-a)*0.5*ξ.+0.5*(a+b)

    int = zeros(ComplexF64,SpaceSize,k+1);
    P = zeros(length(ξ),length(ξ));

    P[:,1] .= 1;

    if k >= 1

        P[:,2] = (2/(b-a)*(t .-a).-1);

        for j = 2:k

          P[:,j+1] = ((2*j-1)*(2/(b-a)*(t .-a).-1).*P[:,j]-(j-1)*P[:,j-1])/j;

        end

    end

    for j = 0:k

      P[:,j+1] *= sqrt((2*j+1)/(b-a));

    end

    for j = 1:(k+1)

        f_j = Assemble_RHS(Nodes,Simplices,Mesh2Space,SpaceSize,N_of_Simplices,Quadrature,f,t[j])

        for i = 1:(k+1)

            int[:,i] += ω[j]*P[j,i]*f_j

        end

    end

    int = 0.5*(b-a)*int

    return reshape(int,(SpaceSize*(k+1),1))

end
