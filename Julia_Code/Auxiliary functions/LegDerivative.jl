function LegDerivative(q::Int64,s::Int64,τ::Float64,side::String)

    if side == "left"

        return (-1)^(s-q)*sqrt((2*s+1)/τ)*(2/τ)^q*derivative_prod(q)*binomial(s+q,s-q)

    elseif side == "right"

        return sqrt((2*s+1)/τ)*(2/τ)^q*derivative_prod(q)*binomial(s+q,s-q)

    end

end
