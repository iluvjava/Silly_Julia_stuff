include("iterative_hessenberg.jl")

# The definition of the Iterative Conjugate Gradient method ====================
#   * only supports real numbers. 
#   * Uses the Iterative Hesseberg that stores the previous 2 orthogonal vectors. 
mutable struct IterativeCGViaArnoldi
    A::Function
    b
    x
    r0
    d  # The previous, last diagonal of the LDL on T
    l_hat  # the last element of the first column of the L^{-1} on L of LDL of T
    l  # the j, j - 1 element of the L matrix from the LDL of T, T is from ih. 
    iter_count::Int64
    ih::IterativeHessenberg
    
    function IterativeCGViaArnoldi(A::Function, b, x0=nothing)
        this = new()
        this.A = A;
        this.b = b;
        this.x = x0 === nothing ? b : x0
        this.r0 = b - A(this.x)
        this.d = 1
        this.l = 1
        this.l_hat = 0
        this.ih = IterativeHessenberg(A, b, x0=x0, max_k=2)
        return this
    end

    function IterativeCGViaArnoldi(A::AbstractArray{T}, b::AbstractArray{T}) where {T<:Real}
        return IterativeCGViaArnoldi((x)-> A*x, b)
    end

end

# Operator Override
function (this::IterativeCGViaArnoldi)()
    ih = this.ih
    if this.iter_count == 0 # The first iteration where conjugate direction is literally r0. 
        ih()
        ih()
        ih()
        H = GetHessenberMatrix(ih)
        T = H[1:2, 1:2]
        display(T)
        return this.x
    end
    q = ih()
    # compute the L matrix for the LDL. 

    this.iter_count += 1
end

