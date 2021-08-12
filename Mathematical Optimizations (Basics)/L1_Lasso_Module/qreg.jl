# Tools for regularized Quantiled regression. 
using COSMO, JuMP, LinearAlgebra
using Statistics
import JuMP
MOI = JuMP.MathOptInterface


function MakeQuantileOptimizationProblem(
    A::Matrix, 
    b::Matrix,
    q::Float64=0.5,
    λ::Float64=0.0
)
    



end








