# This is the original flavor of the Lanczos Algorithm which generates 
# a Tridigonal matrix from a linear operator. 
# * This just performs the Iterative Lanczos, which, will not test to see
# * the problem with matrix, whether is symmetric or not. 

mutable struct IterativeLanczos
    A::Function
    q1
    betas::Dict
    alphas::Dict
    k::Int64
    q_store::Int64
    Q::Dict
    function IterativeLanczos(A::Function, q1, q_store::Int64=typemax(Int64))
        @assert q_store <= 1 "storage must be at least 2. "
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
function (this::IterativeLanczos)()
    
    if this.k == 1
        qNew = this.A(this.Q[this.k])
    else
        qNew = this.A(this.Q[this.k]) - this.betas[this.k - 1]*this.Q[this.k - 1]
    end
    alphaNew = dot(this.Q[this.k], qNew)
    qNew -= alphaNew*this.Q[this.k]
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


function GetQMatrix(this::IterativeLanczos)
    return hcat(values(this.Q)...)
end

function GetTMatrx(this::IterativeLanczos)

    return #TODO: Implement
end

function GetPreviousOrtho(this::IterativeLanczos, j::Int64=0)
    return this.Q[this.k - j]
end