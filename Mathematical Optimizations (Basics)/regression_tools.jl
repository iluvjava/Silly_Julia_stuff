# Using convex solver for Lasso, ridge, and elastic net. 
raw"""
model: y = Xα
y: labels 
X: row data matrix
α: weights
Loss: ‖y - Xα‖^2 + ‖α‖_1

"""
using Convex, SCS, Statistics

function VanderMonde(x::Vector{T}, deg::Int64) where {T<: Number} ::Matrix{T}
    

end


function LassoRidgeElasticNet(
    X::Matrix{T}, 
    y::Matrix{T}, 
    λ::Float64=0.0, # lasso regularizer 
    γ::Float64=0.0, # ridge regularizer
    verbose::Bool=true
) where {T <: Number}
    """
        LassoRidgeElasticNet, it doesn't support biases it will only optimized 
        on weights, sorry about that. 

        Returns: 
            β, the weights, and problem, the optimization object. 

    """
    m, n = size(X)
    @assert size(y, 1) == m && size(y, 2) == 1 "the label and"* 
    " the data matrix doesn't have the maching dimension X:" * 
    string(m, "×",  n) * string(";y: ", size(y))
    @assert γ ≤ 1 && γ ≥ 0 "Regularization γ should be between (0, 1)"
    @assert λ ≤ 1 && λ ≥ 0 "Regularization λ should be between (0, 1)"

    Q = X'*X
    b = y'*X
    β = Variable(n) # the weights 
    
    Loss  = quadform(β, Q)          # β^T*Q*β
    Loss += dot(b, β)               # β^T*b
    Loss += λ*norm(β, 1)            # λ ‖β‖_1
    Loss += γ*sumsquares(β)         # γ ‖β‖_2^2
    
    problem = minimize(Loss)
    solve!(problem, ()-> SCS.Optimizer(verbose=verbose))
    problem.status == SCS.MathOptInterface.OPTIMAL ? vec(evaluate(β)) : β = nothing
    
    return β, problem

end


mutable struct StandardizeRegression
    """
        Standardize so that all features has zero mean, and all labels as zero mean, this is made 
        for the lasso ridge regression function above to take into account of models involved biases. 

    """
    
    X::Matrix # Feature matrix
    y::Matrix # label vector 
    μ::Matrix # feature mean
    ϕ::Number # label mean
    Z::Matrix # Starndardized Matrix
    l::Matrix # Zero mean label 

    function StandardizeIt(X::Matrix, y::Matrix)
        m, _ = size(X)
        @assert size(y, 1) == m && size(y, 2) == 1 "The rows of X should match of the columns of y, but the size of"*
        string("X, Y is: ", size(X), " ", size(y))

        μ = mean(X, dims=1)::Matrix
        ϕ = mean(y)
        Z = X .- μ
        l = y .- ϕ
        new(X, y, μ, ϕ, Z, l)
    end

end

function GetStandardizedData(stdreg::StandardizeRegression)
    return stdreg.X, stdreg.l
end

function Train(stdreg::StandardizeRegression)
    
end

function PredictWithWeights(stadreg::StandardizeRegression, β::Matrix)

    

end

