function LegPols(k::Int64,In::Array{Float64,1},t::Array{Float64,1})

    P = zeros(length(t),k+1);

    a = In[1]; τ = abs(In[2]-a);

    P[:,1] .= 1/sqrt(τ); 

    if k > 0
    
        P[:,2] = sqrt(3/τ)*(2*(t .-a)/τ .-1);

        for j = 3:(k+1)
            P[:,j] = sqrt((2*j-1)/τ)/(j-1)*((2*j-3)*sqrt(τ/(2*j-3))*(2*(t .-a)/τ .-1).*P[:,j-1].-sqrt(τ/(2*j-5))*(j-2)*P[:,j-2]);
        end
    end

    return P
end
