using COSMO, JuMP, LinearAlgebra
using Statistics

function MakeLassoOptimizationProblem(A, y, λ)
    m, n = size(A)
    @assert size(y, 1) == m && size(y, 2) == 1 "the label and"* 
    " the data matrix doesn't have the maching dimension X:" * 
    string(m, "×",  n) * string(";y: ", size(y))
    @assert λ ≤ 1 && λ ≥ 0 "Regularization λ should be between (0, 1)"

    model = Model(with_optimizer(COSMO.Optimizer))
    @variable(model, x[1:n])
    setvalue.(x, A\y)
    @variable(model, η[1:n])
    @constraint(model, -η .<= x)
    @constraint(model, x .<= η)
    @objective(model, Min, λ*sum(η) + sum((A*x - y).^2)) 

    return model

end


mutable struct MyLasso
    ## It's just a collection of data. 

    A::Matrix # Feature matrix
    y::Matrix # label vector 
    μ::Matrix # feature mean
    u::Number # label mean
    Z::Matrix # Starndardized Matrix
    l::Matrix # Zero mean label 
    
    OptModel::Model
    λ::Float64

    function MyLasso(A::Matrix, y::Matrix, λ::Float64=0.0)
        m, _ = size(A)
        @assert size(y, 1) == m && size(y, 2) == 1 ""*
        "The rows of X should match of the columns of y, but the size of"*
        string("X, Y is: ", size(A), " ", size(y))

        μ = mean(A, dims=1)::Matrix
        u = mean(y)
        Z = A .- μ
        l = y .- u
        OptModel = MakeLassoOptimizationProblem(Z, l, λ)
        new(A, y, μ, u, Z, l, OptModel, λ)
    end


end


function LassoAnalyze(this::MyLasso)
    


end


function Changeλ(this::MyLasso, λ)
    """
        Change the Lasso regularizer of the current model 
    """
    model = this.OptModel
    x = model[:x]
    η = model[:η]
    y = this.l
    A = this.Z
    @objective(model, Min, λ*sum(η) + sum((A*x - y).^2)) 
    
end


function SolveForx(this::MyLasso)
    """
        Solve for the weights of the current model 
    """
    optimize!(this.OptModel)
    return value.(this.OptModel[:x])
end

function Getαβ(this::MyLasso)

end