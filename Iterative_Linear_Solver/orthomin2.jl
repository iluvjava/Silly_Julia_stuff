### Implementation of the Orthomin algorithm. It's just conjugate gradient 
### But applied onto the AA^H. 
### This time we use a type so we can interact with all the results. 

using LinearAlgebra
using Logging
mutable struct OrthoMin2
    A::Function
    b
    x0
    p::Dict
    x::Dict 
    r::Dict
    itrCount
    function OrthoMin2(A, b, x0)
        new(A, b, x0, Dict(), Dict(), Dict(), 0)
    end
end

function (this::OrthoMin2)(j::Int64=1)
    @assert j >= 1
    x0 = this.x0
    A = this.A
    b = this.b
    r, x, p = this.r, this.x, this.p
    if this.itrCount == 0
        r[0] = b - A(x0)
        p[0] = r[0]
        x[0] = x0
    end

    for KK ∈ this.itrCount + 1:this.itrCount + j
        Ap = A(p[KK - 1])
        ApAp = Ap⋅Ap
        α = (r[KK - 1]⋅Ap)/(ApAp)
        if α ≈ 0 
            @warn "zero in field of value of matrix A. "
            throw(ErrorException("Compute Error"))
        end
        x[KK] = x[KK - 1] + α*p[KK - 1]
        r[KK] = r[KK - 1] - α*Ap
        β = (A(r[KK])⋅Ap)/(ApAp)
        p[KK] = r[KK] - β*p[KK-1]
    end
    
    Result = x[this.itrCount + j]
    this.itrCount = this.itrCount + j
    return Result
end

function GetResidualNorms(this::OrthoMin2)
    return values(this.r)|>collect.|>x->x⋅x.|>abs.|>sqrt
end

function Test1(N, j)
    # Λ = diagm(rand(N)) .+ 1
    # Q, _ = qr(rand(N, N))
    # A = Q*Λ*transpose(Q)
    A = rand(N, N) + I*(N/10)
    b = rand(N)
    x0 = rand(N)
    this = OrthoMin2(x ->A*x, b, x0)
    this(j)
    println("List of norms for the residual during the iterations")
    for item ∈ this |> GetResidualNorms
        println(item)
    end
end

Test1(40, 30)