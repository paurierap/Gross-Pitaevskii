function ismult(i::Int64,j::Int64)
    
    # Auxiliar function to compute the scalar product P'j*Pi:

    if j < 0 || i > j
        return false
    elseif i == j
        return true
    else
        ismult(i,j-2)
    end   
end
