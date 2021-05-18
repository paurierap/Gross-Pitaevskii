function derivative_prod(x::Int64)

    # Kind of a double factorial function. Needed to calculate the derivative of a Legendre function.

      y = 1;
      for k = 1:2:(2*x-1)
          y *= k;
      end
      return y
end
