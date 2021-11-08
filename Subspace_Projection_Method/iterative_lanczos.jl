# ==============================================================================
# Computing the Lanczos Algorithm
# It will keep the tridiagonal matrix and its LDL decomposition while running. 
# IT will make things easier for another flavor of conjugate Gradient
# * It will store LDL = t
# * By Default it will store all the Q vectors, needs to specify not to store 
#   * them. 
# ==============================================================================

include("iterative_hessenberg.jl")

mutable struct IterativeLanczos
    A
    q
    D              # diagonals of the D matrix for LDL 
    alphas         # diagonals of the T matrix
    betas          # upper/lower diagonals of the T matrix
    L              # Lower diagonals of the L in LDL.T
    Linv           # the first column of L Inverse. 
    store_Q::Int64 # Store all the Orthogonal Vectors
    Q              # Previous 2 orthogonal vectors
    itr::Int64
    function IterativeLanczos(A::Function, q0; store_Q=typemax(Int64))
        if store_Q <= 1
            error("You can do lancozos by storing $(store_Q) number of previous orthogonal vectors")
        end
        this = new()
        this.A = A
        this.q = q0
        this.store_Q = store_Q
        this.Q = Vector{typeof(q0)}()
        push!(this.Q, q0/norm(q0))
        this.alphas = typeof(q0)()
        this.betas = typeof(q0)()
        this.L = typeof(q0)()
        this.D = typeof(q0)()
        this.Linv = typeof(q0)()
        push!(this.Linv, 1)
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
        β = norm(Aq)           # Get the beta here, I assure you it's correct. 
        qNew = Aq/β
        push!(this.Q, qNew)    # Fill in the orthogonal vectors
        push!(this.alphas, α)  # filling in the diagonal for matrix T
        push!(this.betas, β)
        push!(this.D, α)       # fill in the diagonal for D
        this.itr += 1
        return β               # Residual r(1)
    end
    q = this.Q[end]
    Aq = this.A(q)
    Aq -= this.betas[end]*this.Q[end - 1]
    α = dot(q, Aq)
    Aq -= α*q
    β = norm(Aq)
    if β ≈ 0
        error("The q vector is too small after orthogonalization. ")
    end
    qNew = Aq/β

    if length(this.Q) <= this.store_Q
        push!(this.Q, qNew)
    else
        popfirst!(this.Q)
        push!(pushthisQ, qNew)
    end
    push!(this.alphas, α)
    push!(this.L, this.betas[end]/this.D[end])
    push!(this.Linv, -this.Linv[end]*this.L[end])
    d = α - this.betas[end]^2/this.D[end]
    if abs(imag(d)) > 1e-9 || real(2) < 0
        error("The matrix Under Lanczos Might not be Positive Semi-Definite, latest diagonal: $(d)")
    end
    push!(this.D, d)
    push!(this.betas, β)
    
    this.itr += 1
    return β
end

function GetTMatrix(this::IterativeLanczos)
    if this.itr == 0
        error("Cant get the tridiagonal matrix before diagonalizing")
    end
    if this.itr == 1
        return this.alphas[1]
    end
    return SymTridiagonal{Float64}(real(this.alphas), real(this.betas[1:end-1]))
end

function GetQMatrix(this::IterativeLanczos)
    if this.itr == 0 || this.itr == 1
        return this.Q[this.itr + 1]    # Just the fist vector, nothing more. 
    end
    return hcat(this.Q[1:end-1]...)
end

function GetLMatrix(this::IterativeLanczos)
    if this.itr == 0 || this.itr == 1
        return 1
    end
    return Bidiagonal{Float64}(fill(1,this.itr), real(this.L), :L)
end 

function GetDMatrix(this::IterativeLanczos)
    if this.itr == 0 
        error("There is no D matrix yet. Must orthogonalize before getting the D matrix")
    end
    if this.itr == 1
        return this.D[1]
    end
    return Diagonal{Float64}(real(this.D))
end

