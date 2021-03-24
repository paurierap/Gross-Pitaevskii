function derivative_prod(x::Int64)
      y = 1;
      for k = 1:2:(2*x-1)
          y *= k;
      end
      return y
end
