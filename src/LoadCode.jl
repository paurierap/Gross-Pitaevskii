# Load external codes:

include("Assemble_S.jl")
include("./Auxiliary functions/Assemble_Element_S.jl")
include("./Auxiliary functions/LagDerivative.jl")
include("Assemble_M.jl")
include("./Auxiliary functions/Assemble_Element_M.jl")
include("Integrate_Γ.jl")
include("Derivative_Γ.jl")
include("./Auxiliary functions/LegPols.jl")
include("./Auxiliary functions/derivative_prod.jl")
include("./Auxiliary functions/LegDerivative.jl")
include("./Auxiliary functions/Gauss_quad.jl")
include("./Auxiliary functions/ismult.jl")
include("Fixed_Point.jl")
include("Mixed_BC.jl")
include("RHS_Integration.jl")
include("mass_integration.jl")
include("energy_integration.jl")
include("./Auxiliary functions/GRnw.jl")
include("./Auxiliary functions/GLnw.jl")
include("Gr.jl")