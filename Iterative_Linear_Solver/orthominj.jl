### Orthogonalizing all previous direction using modified gram-schidmtz. 
### And then use them to reduce the residual vector under the 2-norm. 
###   * Store the MGS subspace for the Krylov Subspace 
###   * Support iterators
### Algorithm Features: 
###   * Proposed as orthomin(j) methods in Greenbaum's book. It compute krylov 
###   * subspace κ and A(κ) on the go and orthogonalizing it on the go. 

using Logging
using LinearAlgebra

mutable struct OrthoMin
    
    A::Function # linear operator, please don't passing mutating linear operator. 
    b
    x0
    itr::Int64
    p::Vector
    s::Vector
    x::Vector
    r::Vector
    storageSize::Int64

    function OrthoMin(A::Function, b, x0; storage_size=10)
        this = new()
        this.A = A
        this.b = b
        this.x0 = x0
        this.itr = 0
        this.p = Vector()
        this.s = Vector()
        this.x = Vector()
        this.r = Vector()
        this.storageSize = storage_size
        return this
    end

    function OrthoMin(A::AbstractMatrix, b, x0; storage_size=10)
        return OrthoMin((x)-> (A*x), b, x0; storage_size=10)
    end
    
end


function (this::OrthoMin)()
    # Does one exact iteration of orthomin j. 
    prj(x, y) = (dot(x,y))/(dot(y, y))
    A, b, x0 = this.A, this.b, this.x0
    r, x, p, s = this.r, this.x, this.p, this.s

    if this.itr == 0  # initialization
        push!(r, b - A(x0))
        push!(x, x0)
        push!(p, r[1])
        push!(s, A(r[1]))
    end
    a = prj(r[end], s[end])
    println(a)
    if a ≈ 0
        @warn "Field of value contains zero."
        throw(ErrorException("Computational Error"))
    end
    if isinf(a)
        @warn "a approaching infinity, s, direction vector is way too small"
        throw(ErrorException("Computational Error"))
    end
    push!(x, x[end] + a*p[end])
    push!(r, r[end] - a*s[end])
    pnew = copy(s[end])
    snew = A(s[end])
    for LL ∈ 1:length(s)
        β = prj(snew, s[LL])
        pnew -= β*p[LL]
        snew -= β*s[LL]
    end
    push!(p, pnew)
    push!(s, snew)
    
    if this.storageSize <= 0 && length(x) > this.storageSize
        popfirst!(x)
        popfirst!(r)
        popfirst!(p)
        popfirst!(s)
    end
    this.itr += 1
    return x[end]

end

function GetAllResidualMeasure(this::OrthoMin)
    return this.r .|> x -> x⋅x .|> abs .|> sqrt
end

function Test1(N)
    A = diagm(rand(N) .- 0.5)
    b = rand(N)
    x0 = (A\b) + randn(N)*1e-1
    Instance = OrthoMin(A, b, x0, storage_size=-1)
    for II ∈ 1:N
        Instance()
    end
    println("List of residual norm is: ")
    Instance|>GetAllResidualMeasure|>display
end

Test1(10)