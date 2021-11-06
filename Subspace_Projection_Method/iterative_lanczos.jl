# ==============================================================================
# Computing the Lanczos Algorithm
# It will keep the tridiagonal matrix and its LDL decomposition while running. 
# IT will make things easier for another flavor of conjugate Gradient
# * It will store LDL = t
# * It won't store the matrix Q
# ==============================================================================

include("iterative_hessenberg.jl")

mutable struct IterativeLanczos
    A
    q
    D
    T
    L
    store_Q::Bool
    Q
    function IterativeLanczos(A::Function, q0; store_Q::Bool=false)
        this = new()
        this.A = A
        this.store_Q = store_Q
        this.Q = vector{typeof(q0)}()
        push!(this.Q, q0)
        return this
    end

    function IterativeLanczos(A::AbstractArray, q0)
        return IterativeLanczos((x) -> A*x, q0)
    end

end