include("iterative_hessenberg.jl")

# The definition of the Iterative Conjugate Gradient method ====================


mutable struct IterativeConjugateGradient
    A::Function
    b
    x0
    iter_count::Int64
    function IterativeConjugateGradient(f::Function, b, x0)
        this = new()
        this.A = f;
        this.b = b;
        this.x0 = x0  
    end

    function IterativeConjugateGradient(A::AbstractArray{T}, b::AbstractArray{T}) where {T<:Real}
        
    end

end


# Operator Override
function (this::ConjugateGradient)()

end

