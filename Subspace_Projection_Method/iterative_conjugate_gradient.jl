include("iterative_hessenberg.jl")

# The definition of the Iterative Conjugate Gradient method ====================
#   * only supports real numbers. 
#   * Uses the Iterative Hesseberg that stores the previous 2 orthogonal vectors. 
mutable struct IterativeConjugateGradient
    A::Function
    b
    x
    r0
    iter_count::Int64
    l1 # the j, j - 1 element of the L matrix from the LDL of T, T is from ih. 
    function IterativeConjugateGradient(f::Function, b, x0=nothing)
        this = new()
        this.A = f;
        this.b = b;
        this.x = x0 === nothing ? b : x0
        this.r0 = b - A(this.x)
    end

    function IterativeConjugateGradient(A::AbstractArray{T}, b::AbstractArray{T}) where {T<:Real}
        IterativeConjugateGradient((x)-> A(x), b)
    end

end

# Operator Override
function (this::IterativeConjugateGradient)()
    ih = this.ih
    if this.iter_count == 0 # The first iteration where conjugate direction is literally r0. 
        q = ih()
        this.x += (ih.H[end][end]/ih.H[end][end - 1])*q
        this.iter_count += 1
        return ih.H[end][end]
    end
    q = ih()
    # compute the L matrix for the LDL. 

    this.iter_count += 1
end

