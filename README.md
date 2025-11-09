# 1-D Gross-Pitaevskii equation solver 

## Summary

- The Finite Element Method is used for space discretization. The solver supports arbitrary high order Lagrange polynomials (typified by the parameter `d`).

- A collocation Method described in [[1]](#1) is employed for time integration (with parameters `k` and `r`). 

## Usage

To run, create a spatial mesh using `Generate_Mesh1D.jl`, which also provides a time discretization for $(0,T]$ given in a `.jld` file. Then run `main.jl`, which solves the GPE with $\beta=-2$, providing a `.jld` file with the approximate solution at $t=T$, along with the values of the discretized energy and mass. 

## References
<a id="1">[1]</a>
S. Becher, G. Matthies, and D. Wenzel (2018). *Variational methods for stable time discretization of first-order differential equations*. doi: https://www.researchgate.net/publication/322821295_Variational_methods_for_stable_time_discretization_of_first-order_differential_equations

