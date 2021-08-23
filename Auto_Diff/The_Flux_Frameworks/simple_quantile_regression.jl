using Flux

mutable struct QuantileRegression

    W::VecOrMat
    b::Union{Vector, Number}
    q::Number

    # the gradient tape
    g::Flux.Zygote.Grads

    

    function QuantileRegression(features_in::Int64, dim_out::Int64; q::Union{Number, Nothing}=nothing)
        W = rand(features_in, dim_out)
        b = rand(dim_out)
        this = new()
        this.W = W
        this.b = b 
        this.q = q === nothing ? 0.5 : q
        return this
    end

end


function Predict(this::QuantileRegression, X::Matrix)
    return X*this.W .+ this.b
end

function Loss(this::QuantileRegression, X::Matrix, y::VecOrMat)
    ŷ = Predict(this, X)
    α = this.q
    l = map(ŷ - y) do d
        if d > 0 
            return α*d
        else
            return (α - 1)*d
        end
    end
    return sum(l)
end

"""
    Get the gradient, give the data X, and labels y. 
"""
function ∇(this::QuantileRegression, X::Matrix, y::VecOrMat)    
    this.g = gradient(()-> Loss(this, X, y), params(this.W, this.b))
    this.X = X
    this.y = y
    return [this.g[p] for p ∈ this.g.params]
end


"""
    Train the model with classic momentum and armijo line search, adaptive 
    reset and stepsize adjustment. 
"""
function Train!(
        this::QuantileRegression,
        X::Matrix, 
        y::VecOrMat; 
        η::Float64=1e-3, 
        α::Float64=0.2
    )

    
end