function GLnw(k::Int64)

    N1 = k+1;

    # Use the Chebyshev-Gauss-Lobatto nodes as the first guess.

    x = cos.(π*(0:k)/k);

    # The Legendre Vandermonde Matrix

    P = zeros(N1,N1);

    # Compute P_(N) using the recursion relation. Compute its first and second derivatives and
    # update x using the Newton-Raphson method.

    xold = 2.;

    while maximum(abs.(x.-xold)) > 1e-12
        xold = x;
        P[:,1] .= 1.;    P[:,2] = x;
        for j = 2:k
            P[:,j+1] = ((2*j-1)*x.*P[:,j]-(j-1)*P[:,j-1])/j;
        end
        x = xold .- (x.*P[:,N1]-P[:,k])./(N1*P[:,N1]);
    end

    w = 2 ./(k*N1*P[:,N1].^2);

    return x,w
end
