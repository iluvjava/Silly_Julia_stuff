include("iterative_lanczos.jl")

# The definition of the Iterative Conjugate Gradient method ====================
#   * only supports real numbers. 
#   * Uses the Iterative Lanczos with LDL partially stored. 
#   
mutable struct IterativeCGViaLanczos
    A::Function
    b      # The RHS of the equation. 
    x      # The solution vector. 
    r0     # The first initial residual. 
    r     # That norm of the newest residual. 
    p      # The last conjugate directions. 
    itr::Int64
    il::IterativeLanczos
    
    function IterativeCGViaLanczos(A::Function, b, x0=nothing)
        this = new()
        this.A = A;
        this.b = b;
        this.x = x0 === nothing ? b : x0
        this.r0 = b - A(this.x)
        this.r = norm(this.r0)
        this.il = IterativeLanczos(A, this.r0, store_Q=2)
        this.p = this.il.Q[1]  # conjugate direction directly come from Q from lanczos
        this.itr = 0
        return this
    end

    function IterativeCGViaLanczos(A::AbstractArray{T}, b::AbstractArray{T}) where {T<:Real}
        return IterativeCGViaLanczos((x)-> A*x, b)
    end

end

# Operator Override
function (this::IterativeCGViaLanczos)()
    il = this.il
    b = this.b
    A = this.A
    if this.itr == 0 # The first iteration where conjugate direction is literally r0. 
        β = il()
        a = 1/il.alphas[end]  # step size
        this.x += a*this.r0
        this.r = this.r0 - a*A(this.r0)  # TODO: Is it the same as the Lanczos vectors? 
        this.itr += 1
        return a*β  # 2 norm of the residual. 
    end
    # q = ih()
    # compute the L matrix for the LDL. 
    

    this.itr += 1
    return nothing
end

A = rand(3,3)
A = A*A'
b = rand(3)
cg = IterativeCGViaLanczos(A, b)
rnorm = cg();
il = cg.il
