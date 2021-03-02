function Fixed_Point(β,b,x_0,A,Nodes::Array{Float64,1},Simplices::Array{Int64,2},Mesh2Space::Array{Int64,1},SpaceSize::Int64,k::Int64,r::Int64,In)

    Δx = ones(SpaceSize)

    x_next = x_0

    it = 0

    τ = abs(In[2]-In[1]) # timestep.

    if r == 1  #cGP(k)

        while (norm(Δx) > 1e-12)

            x_n = copy(x_next)

            C = copy(A)

            for i = 0:(k-1)

                for j = 0:k

                    A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] -= β*Integrate_M_tau(x_n,i,j,k,In,Nodes,Simplices,Mesh2Space,SpaceSize)
                    
                end

            end

            x_next = A\b;

            it += 1

            Δx = x_next-x_n

            A = C

        end

    else  #G_r(k)

        D0 = floor(Int,0.5*(r-1))
        D1 = floor(Int,0.5*r)

        while (norm(Δx) > 1e-9)

            x_n = copy(x_next)

            C = copy(A)

            # In the local problem, the matrix M_tau needs to be integrated:

            for i = 0:(k-r)

                for j = 0:k

                    A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] -= β*Integrate_M_tau(x_n,i,j,k,In,Nodes,Simplices,Mesh2Space,SpaceSize)

                end

            end

            # Finally, continuity conditions on the derivatives depending on D0 and D1, where M_tau is NOT integrated:

            for q = 1:D0

                i = k-r+1+q

                side = "left"

                for j = 0:k

                    for α = 0:(q-1)

                        A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] -= β*binomial(q-1,α)*LegDerivative(α,j,τ,side)*Derivative_M_tau(q-1-α,k,τ,x_n,Nodes,Simplices,Mesh2Space,SpaceSize,side)

                    end

                end

            end

            for q = 1:D1

                i = k-r+1+q+D0

                side = "right"

                for j = 0:k

                    for α = 0:(q-1)

                        A[i*SpaceSize+1:(i+1)*SpaceSize,j*SpaceSize+1:(j+1)*SpaceSize] -= β*binomial(q-1,α)*LegDerivative(α,j,τ,side)*Derivative_M_tau(q-1-α,k,τ,x_n,Nodes,Simplices,Mesh2Space,SpaceSize,side)

                    end

                end

            end

            x_next = A\b; # Solve to get the Legendre coefficients to later check convergence.

            it += 1

            Δx = x_next-x_n

            A = C

        end

        # Usually converges in 6-8 iterations!

    end

    return x_next

end
