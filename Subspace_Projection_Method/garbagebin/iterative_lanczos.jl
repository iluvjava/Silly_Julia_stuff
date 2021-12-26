# This is the original flavor of the Lanczos Algorithm which generates 
# a Tridigonal matrix from a linear operator. 
#   * This just performs the Iterative Lanczos, which, will not test to see
#   * the problem with matrix, whether is symmetric or not. 

mutable struct IterativeLanczos
    A::Function
    q1
    betas::Dict
    alphas::Dict
    k::Int64       # when k=1, the algorithm haven't started yet. 
    q_store::Int64
    Q::Dict
    function IterativeLanczos(A::Function, q1, q_store::Int64=typemax(Int64))
        @assert q_store > 1 "storage must be at least 2. "
        this = new()
        this.q1 = copy(q1)
        this.A = A
        this.betas = Dict()
        this.alphas = Dict()
        this.k = 1
        this.q_store = q_store
        this.Q = Dict{Int64, typeof(q1)}()
        this.Q[1] = q1/norm(q1)
        return this
    end

    function IterativeLanczos(A::AbstractArray, q1, q_store::Int64=typemax(Int64))
        return IterativeLanczos((x) -> A*x, q1, q_store)
    end

end


# Override the () operator. 
"""
    It performs one iteration of the Lanczos Orthogonalization. 
"""
function (this::IterativeLanczos)()
    q = this.Q[this.k]
    Aq = this.A(q)
    if this.k == 1
        qNew = Aq
    else
        qPre = this.Q[this.k - 1]
        qNew = Aq - this.betas[this.k - 1]*qPre
    end
    alphaNew = dot(q, qNew)
    qNew -= alphaNew*q
    betaNew = norm(qNew)
    qNew /= betaNew

    if length(this.Q) > this.q_store
        delete!(this.Q, this.k + 1 - length(this.Q))
    end

    this.Q[this.k + 1] = qNew
    this.alphas[this.k] = alphaNew
    this.betas[this.k] = betaNew
    this.k += 1
    return betaNew
end

"""
    It calls the () operator 3 j times, performing j iterations. 
"""
function (this::IterativeLanczos)(j::Int64)
    betas = Vector()
    for _ in 1: j
        push!(betas, this())
    end
    return betas
end


"""
    Get the Q matrix from the Lanczos Iterations, the number of columnos of the 
    matrix depends on this.q_store.
    * k by k matrix returned. 

"""
function GetQMatrix(this::IterativeLanczos)
    Q = this.Q
    qs = Vector{typeof(Q[1])}()
    for Idx in 1: this.k
        push!(qs, this.Q[Idx])
    end
    return hcat(qs...)
end


# """
#     Get the H matrix, the Hessenberg form, it will be k by k - 1.
    
# """
# function GetHMatrix(this::IterativeLanczos)
#     x, y, values = Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
#     for Idx in 1:this.k

#     end
# end


"""
    Return the Tridiagonal matrix, take note thta, the matrix only supports 
    floating point, and it's returned as a sparse matrix. 
    * A (k - 1) by (k - 1) matrix is returned. 

"""
function GetTMatrix(this::IterativeLanczos)
    if this.k == 1
        return nothing
    end
    if this.k == 2
        return this.alphas[1]
    end
    return SymTridiagonal{Float64}(
        [real(this.alphas[Idx]) for Idx in 1: this.k - 1],
        [real(this.betas[Idx]) for Idx in 1:this.k - 2]
    )
end

"""
    Get the previous 3 Orthogonal vector from the Lanczos Algorithm. 
    * If we haven' generated 3 orthogonal vector 3, then 
    this will just return the first 2, or 1 orthgonal matrix. 
"""
function GetPrevious3OrthogonalVec(this::IterativeLanczos)
    if this.k == 1
        return this.Q[1]
    end
    if this.k == 2
        return hcat(this.Q[1], this.Q[2])
    end
    return hcat(this.Q[this.k], this.Q[this.k - 1], this.Q[this.k - 2])
end