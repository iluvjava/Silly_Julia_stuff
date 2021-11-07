include("iterative_hessenberg.jl")  # import stuff that this import. 

mutable struct IterativeCGOriginal
    r
    rnew
    d
    A::Function
    x
    b
    itr
    function IterativeCGOriginal(A::Function, b, x0=nothing)
        this = new()
        this.A = A
        this.x = x0 === nothing ? b : x0
        this.r = b - A(this.x)
        this.rnew = similar(this.r)
        this.d = this.r
        this.itr = 0
        return this
    end

    function IterativeCGOriginal(A::AbstractArray, b::AbstractArray)
        return IterativeCGOriginal((x)->A*x, b)
    end
    
end

function (this::IterativeCGOriginal)()
    r = this.r
    A = this.A
    d = this.d
    Ad = A(d)

    α = dot(r, r)/dot(d, Ad)
    if α < 0 
        error("CG got a non-definite matrix")
    end
    this.x += α*d
    this.rnew = r - α*Ad                # update rnew 
    β = dot(this.rnew, this.rnew)/dot(r, r)
    # @assert abs(dot(rnew + β*d, Ad)) < 1e-8 "Not conjugate"
    this.d = this.rnew + β*d
    this.r = this.rnew                      # Override
    this.itr += 1 
    return r
end

