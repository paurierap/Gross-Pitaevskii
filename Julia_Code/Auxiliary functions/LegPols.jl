function LegPols(k::Int64,In::Array{Float64,1},t::Array{Float64,1})

    P = zeros(length(t),k+1);

    a = In[1]; b = In[2];

    τ = abs(b-a);

    P[:,1] .= 1.; P[:,2] = (2*(t .-a)/τ .-1);

    for j = 3:(k+1)
        P[:,j] = 1/j*((2*j-1)*(2*(t .-a)/τ .-1).*P[:,j-1].-(j-1)*P[:,j-2]);
    end

    for j = 1:k+1
        P[:,j] *= sqrt((2*j-1)/τ);
    end
    
    return P
end
