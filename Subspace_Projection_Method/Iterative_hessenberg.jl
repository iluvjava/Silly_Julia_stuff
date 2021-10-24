using LinearAlgebra

# ==============================================================================
# ITERATIVE HESSENBERG
# ==============================================================================
mutable struct IterativeHessenberg
    A::Function
    b
    H::Vector  # vector of vector
    Q::Vector  # vector of orthogonal vector
    r0
    x0
    
    function IterativeHessenberg(A::Function, b, x0)
        this = new()
        this.A = A
        this.b = b
        this.x0 = x0
        this.r0 = b - A(x0)
        this.H = Vector{Vector{Number}}()
        this.Q = Vector{typeof(b)}()
        push!(this.Q, this.r0./sqrt(dot(this.r0, this.r0)))
        return this
    end

    function IterativeHessenberg(A::Matrix, b::VecOrMat)
        this = IterativeHessenberg((x) -> A*x, b, b)
        return this
    end

end

## Functional Operator Override. 
function (this::IterativeHessenberg)()
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
    return Q[end], this.H[end]
end

# Access elements and objects from this algorithms. 

function GetHessenberMatrix(this::IterativeHessenberg)
    n = length(this.H)
    m = n + 1
    H = Matrix{typeof(this.H[1][1])}(undef, m, n)
    for IdxJ ∈ 1:n
        for IdxI ∈ 1:length(this.H[IdxJ])
            H[IdxI, IdxJ] = this.H[IdxJ][IdxI]
        end
        for IdxI ∈ length(this.H[IdxJ])+1:m 
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


# ==============================================================================


# A brief test. 
A = rand(3,3)
ib = IterativeHessenberg(A, rand(3,1))
ib()
ib()
ib()
Q = GetOrthogonalMatrix(ib)
println("The Q matrix is: ")
display(Q)
H = GetHessenberMatrix(ib)
println("The H matrix is: ")
display(H)
println("Testing the Recurrence: AQ_{k-1} = Q_{k}H") 
println("AQ_{k-1} is: ")
display(A*Q[:, 1:end-1])
println("Q_{k}H is:")
display(Q*H)
