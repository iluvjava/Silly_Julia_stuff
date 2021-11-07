include("iterative_lanczos.jl")

# The definition of the Iterative Conjugate Gradient method ====================
#   * only supports real numbers. 
#   * Uses the Iterative Hesseberg that stores the previous 2 orthogonal vectors. 
mutable struct IterativeCGViaArnoldi
    A::Function
    b
    x
    r0
    p      # The last conjugate directions
    iter_count::Int64
    il::IterativeLanczos
    
    function IterativeCGViaArnoldi(A::Function, b, x0=nothing)
        this = new()
        this.A = A;
        this.b = b;
        this.x = x0 === nothing ? b : x0
        this.r0 = b - A(this.x)
        this.p = r0/norm(r0)  # matches with lancozos q vector. 
        this.il = IterativeLanczos(A, this.r0, store_Q=2)
        return this
    end

    function IterativeCGViaArnoldi(A::AbstractArray{T}, b::AbstractArray{T}) where {T<:Real}
        return IterativeCGViaArnoldi((x)-> A*x, b)
    end

end

# Operator Override
function (this::IterativeCGViaArnoldi)()
    il = this.il
    if this.iter_count == 0 # The first iteration where conjugate direction is literally r0. 
        il()
        a = dot(this.p, this.p)/il.betas[1]  # step size
        this.x -= a*norm(this.r0)*this.p
        return this.alphas[1]  # 2 norm of the residual. 
    end
    # q = ih()
    # compute the L matrix for the LDL. 

    this.iter_count += 1
end

