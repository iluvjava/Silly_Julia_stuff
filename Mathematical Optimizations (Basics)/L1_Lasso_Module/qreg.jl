# Tools for regularized Quantiled regression. 
using COSMO, JuMP, LinearAlgebra
using Statistics
MOI = MathOptInterface


function MakeQuantileOptimizationProblem(
    A::Matrix, 
    b::Union{Matrix, Vector},
    q::Float64=0.5,
    λ::Float64=0.0, 
    silent::Bool=true
)
    m, n = size(A)
    if isa(x, Vector)
        x = reshape(x, (length(x), 1))
    end
    @assert size(x, 1) == size(A, 2) "The number of rows in A should be the"*
    "The same as the height of the vector x, but, A is $(size(A)) and b is $(size(b))"
    @assert q >= 0 && q <= 1 "The quantile should be between 0 and 1"*
    "don't use precentage please"

    # pass constructor to JuMP interface
    model = Model(with_optimizer(COSMO.Optimizer))  
    set_optimizer_attribute(model, MOI.Silent(), silent) 
    @variable(model, x[1:n])
    @variable(model, η[1:n])
    @variable(model, ϵ[1:n])
    
    @constraints(model, begin
        x .<= η
        η .>= x
        A*x - b .>= q*ϵ
        -(1 - q)*ϵ .<= A*x - b
    end)

    @objective(model, Min, λ*sum(η) + sum(ϵ))

    return model
end








