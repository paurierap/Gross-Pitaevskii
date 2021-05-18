function LegDerivative(q::Int64,s::Int64,τ::Float64,side::String)

    # Calculate the q-th derivative of the s+1-th Legendre Polynomial at either side.

    if side === "left"
        return (-1)^(s-q)*sqrt((2*s+1)/τ)*(2/τ)^q*derivative_prod(q)*binomial(s+q,s-q)
    end

    return sqrt((2*s+1)/τ)*(2/τ)^q*derivative_prod(q)*binomial(s+q,s-q)
    
end
