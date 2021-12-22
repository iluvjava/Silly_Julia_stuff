using LinearAlgebra

# The definition of the Iterative Conjugate Gradient method ====================
#   * only supports real numbers. 
#   * Uses the Iterative Lanczos with LDL partially stored. 
#   * It doen't store the residual because of a sign problem between lanczos and
#   original cg. 
#   * not sure how to get the one norm of the residual. 

mutable struct IterativeCGViaLanczos
    A::Function
    b      # The RHS of the equation. 
    x      # The solution vector, when it's all the guesses including the initial guess. 
    r0     # The first initial residual. 
    r0norm
    r      # That norm of the newest residual. 
    p      # The last conjugate directions. 
    itr::Int64
    il::IterativeLanczos
    
    function IterativeCGViaLanczos(A::Function, b, x0=nothing)
        this = new()
        this.A = A;
        this.b = b;
        this.x = x0 === nothing ? b .+ 0.1 : x0   # perturb initial guess to handle eigenvalue of 1 of the matrix A. 
        this.r0 = b - A(this.x)
        this.r = norm(this.r0)
        this.r0norm = this.r
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
    if this.itr == 0 # The first iteration where conjugate direction is literally r0. 
        β = il()
        a = 1/il.alphas[end]
        this.x += a*this.r0
        this.r = abs(this.r*β*a)
        this.itr += 1
        return this.r 
    end
    β = il()
    a = (1/il.D[end])*il.Linv[end]*this.r0norm
    this.p = il.Q[end - 1] - il.L[end]*this.p
    this.x += a*this.p
    this.r = abs(this.r*β*a)
    this.itr += 1
    return this.r
end



