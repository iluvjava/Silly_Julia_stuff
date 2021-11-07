# ==============================================================================
# Computing the Lanczos Algorithm
# It will keep the tridiagonal matrix and its LDL decomposition while running. 
# IT will make things easier for another flavor of conjugate Gradient
# * It will store LDL = t
# * It won't store the matrix Q with more than 2 vectors init. 
# ==============================================================================

include("iterative_hessenberg.jl")

mutable struct IterativeLanczos
    A
    q
    D              # diagonals of the D matrix for LDL 
    alphas         # diagonals of the T matrix
    betas          # upper/lower diagonals of the T matrix
    L              # Lower diagonals of the L in LDL.T
    store_Q::Bool  # Store all the Orthogonal Vectors
    Q              # Previous 2 orthogonal vectors
    itr::Int64
    function IterativeLanczos(A::Function, q0; store_Q=false)
        this = new()
        this.A = A
        this.q = q0
        this.store_Q = store_Q
        this.Q = Vector{typeof(q0)}()
        push!(this.Q, q0/norm(q0))
        this.alphas = Vector{Real}()
        this.betas = Vector{Real}()
        this.L = Vector{Real}()
        this.D = Vector{Real}()
        this.itr = 0
        return this
    end

    function IterativeLanczos(A::AbstractArray, q0)
        return IterativeLanczos((x) -> A*x, q0)
    end

end


function (this::IterativeLanczos)() 
    if this.itr == 0
        q = this.Q[end] 
        Aq = this.A(q)
        α = dot(q, Aq)
        Aq -= α*q              # Remove the projection onto q. 
        β = norm(Aq)
        qNew = Aq/β
        
        push!(this.Q, qNew)    # Fill in the orthogonal vectors
        push!(this.alphas, α)  # filling in the diagonal for matrix T
        push!(this.betas, β)
        push!(this.D, α)       # fill in the diagonal for D

        this.itr += 1
        return β               # Residual r(1)
    end

    this.itr += 1
    return 
end

function GetTridiagonalMatrix(this::IterativeLanczos)
    if this.itr == 0
        error("Cant get the tridiagonal matrix before diagonalizing")
    end
    return SymTridiagonal(this.alphas, this.betas)
end