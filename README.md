# Gross-Pitaevskii 

1-D Gross-Pitaevskii equation solver: FEM discretization in space (d-order polynomials) and collocation method for time integration (with parameters k,r). To run, create a spatial mesh using ```Generate_Mesh1D.jl```, which also provides a time discretization from ![formula](https://render.githubusercontent.com/render/math?math=\(0,T\]). Then run ```main.jl```, which solves the GPE with ![formula](https://render.githubusercontent.com/render/math?math=\beta=-2).
