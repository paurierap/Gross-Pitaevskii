function GRnw(k::Int64)

    N1 = k+1;

    # Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes.

    x = -cos.(2*π*(0:k)/(2*k+1));

    # The Legendre Vandermonde Matrix

    P = zeros(N1,N1+1);

    # Compute P_(N) using the recursion relation. Compute its first and second derivatives and
    # update x using the Newton-Raphson method.

    xold = 2;

    # Free abscissae

    free = 2:N1;

    while maximum(abs.(x.-xold)) > 1e-12

        xold = x;

        P[1,:] = (-1).^(0:N1);

        P[free,1] .= 1;    P[free,2] = x[free];

        for j = 2:N1

          P[free,j+1]=((2*j-1)*x[free].*P[free,j]-(j-1)*P[free,j-1])/j;

        end

        x = [-1;xold[free]-((1 .-xold[free])/N1).*(P[free,N1]+P[free,N1+1])./(P[free,N1]-P[free,N1+1])];

    end

    # The Legendre-Gauss-Radau Vandermonde

    P = P[1:N1,1:N1]

    # Compute the weights

    w = zeros(N1,1);
    w[1] = 2/N1^2;
    w[free]=(1 .-x[free])./(N1*P[free,N1]).^2;

    return x,w

end
