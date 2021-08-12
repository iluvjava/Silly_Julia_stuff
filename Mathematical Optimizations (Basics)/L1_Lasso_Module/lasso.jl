using COSMO, JuMP, LinearAlgebra
using Statistics
import JuMP
MOI = JuMP.MathOptInterface

function MakeLassoOptimizationProblem(A::Matrix, y::Matrix, λ::Float64)
    """
        Phrase the quadratic programming problem for Lasso regularization 
        problem. 

    """
    m, n = size(A)
    @assert size(y, 1) == m && size(y, 2) == 1 "the label and"* 
    " the data matrix doesn't have the maching dimension X:" * 
    string(m, "×",  n) * string(";y: ", size(y))
    @assert λ ≤ 1 && λ ≥ 0 "Regularization λ should be between (0, 1)"

    model = Model(with_optimizer(COSMO.Optimizer))

    set_optimizer_attribute(model, MOI.Silent(), true)
    set_optimizer_attribute(model, "max_iter", 50000)

    @variable(model, x[1:n])
    setvalue.(x, A\y)
    @variable(model, η[1:n])
    @constraint(model, -η .<= x)
    @constraint(model, x .<= η)
    @objective(model, Min, λ*sum(η) + sum((A*x - y).^2)) 

    return model

end


mutable struct LassoSCOP
    ## It's just a collection of data. 

    A::Matrix # Feature matrix
    y::Matrix # label vector 
    μ::Matrix # feature mean
    u::Number # label mean
    Z::Matrix # Starndardized Matrix
    l::Matrix # Zero mean label 
    
    OptModel::Model
    λ::Float64

    function LassoSCOP(A::Matrix, y::Matrix, λ::Float64=0.0)
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


function LassoPath(this::LassoSCOP)
    """
        Analyze the Lasso Problem. 

    """

    u = this.u
    A = this.Z
    y = this.l
    Results = Vector{Vector}()
    
    function λMax(A, y)
        ToMax = A'*(y .- u)
        ToMax *= 2
        return maximum(abs.(ToMax))
    end

    λ = λMax(A, y)
    Changeλ(this, λ)
    x = SolveForx(this)
    push!(Results, x)
    MaxItr = 100

    while λ >= 1e-16 && MaxItr >= 0
        λ /= 2
        MaxItr -= 1
        setvalue.(this.OptModel[:x], value.(x)) 
        Changeλ(this, λ)
        x = SolveForx(this)
        push!(Results, value.(x))
    end

    ResultsMatrix = zeros(length(x), length(Results))
    for II ∈ 1:length(Results)
        ResultsMatrix[:, II] = Results[II]
    end

    return ResultsMatrix
end


function Changeλ(this::LassoSCOP, λ)
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


function SolveForx(this::LassoSCOP)
    """
        Solve for the weights of the current model, given the current configuration of the model.
    """
    optimize!(this.OptModel)
    
    TernimationStatus = termination_status(this.OptModel)
    @assert TernimationStatus == MOI.OPTIMAL "Terminated with non-optimal value when solving for x. "*
    string("The status is: ", TernimationStatus)
    return value.(this.OptModel[:x])
end


function Getαβ(this::LassoSCOP, lambda::Float64)
    """
        Get the weights and biases, for the original model (The model trained is normalized), for a 
        particular regularization value. 
    """

    
    
end


