using LinearAlgebra

# ==============================================================================
# ITERATIVE HESSENBERG
#   * Type of A, and type of b has to be the same. For example, if A is complex, 
#   then b vector/matrix has to be complex too. 

mutable struct IterativeHessenberg
    
    A::Function
    b
    H::Vector  # vector of vector
    Q::Vector  # vector of orthogonal vector
    r0
    x0
    maxK::Int64
    itr_count::Int64
    
    function IterativeHessenberg(A::Function, b; x0=nothing, max_k=typemax(Int64))
        this = new()
        this.A = A
        this.b = b
        this.x0 = x0 === nothing ? b : x0
        this.r0 = b - A(this.x0)
        this.H = Vector{Vector{Number}}()
        this.Q = Vector{typeof(b)}()
        push!(this.Q, this.r0./sqrt(dot(this.r0, this.r0)))
        this.maxK = max_k
        this.itr_count = 0
        return this
    end

    function IterativeHessenberg(
            A::Matrix, 
            b::VecOrMat; 
            x0=nothing, 
            max_k=typemax(Int64)
        )
        this = IterativeHessenberg((x) -> A*x, b; x0=x0, max_k=max_k)
        return this
    end

end

## Functional Operator Override. 
function (this::IterativeHessenberg)()
    if this.itr_count == 0
        this.itr_count += 1
        return this.r0, sqrt(dot(this.r0, this.r0))
    end
    A = this.A
    Q = this.Q
    u = A(Q[end])
    h = Vector()
    for q ∈ Q
        push!(h, dot(q, u))
        u .-= h[end].*q
    end
    push!(h, sqrt(dot(u, u)))
    push!(this.H, h)
    push!(Q, u./h[end])
    
    if length(Q) > this.maxK
        popfirst!(Q)
    end
    this.itr_count += 1
    return Q[end], this.H[end][end]
end

# Access elements and objects from this algorithms. 

function GetHessenberMatrix(this::IterativeHessenberg)
    n = length(this.H)
    m = n + 1
    H = Matrix{typeof(this.H[1][1])}(undef, m, n)
    for IdxJ ∈ 1:n
        for IdxI ∈ (IdxJ + 1):-1:1
            OffsetFromDiag = IdxJ + 1 - IdxI
            if OffsetFromDiag < length(this.H[IdxJ])
                H[IdxI, IdxJ] = this.H[IdxJ][end - OffsetFromDiag]
            else
                H[IdxI, IdxJ] = 0
            end
        end
        for IdxI ∈ (IdxJ + 2):m  # zeros below everything. 
            H[IdxI, IdxJ] = 0
        end
    end
    
    return H
end

function GetOrthogonalMatrix(this::IterativeHessenberg)
    Q = Matrix{typeof(this.Q[1][1])}(undef, length(this.Q[1]), length(this.Q))
    for IdxJ ∈ 1:size(Q, 2)
        Q[:, IdxJ] = this.Q[IdxJ]
    end
    return Q
end


