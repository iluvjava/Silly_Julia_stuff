include("iterative_hessenberg.jl")

# The instane of conjugate gradient method. 

mutable struct ConjugateGradient
    ih:: IterativeHessenberg
    A::Function
    residual
    x
    b

    function ConjugateGradient(A::Function, b; x0=nothing)
        this = new()
        this.ih = IterativeHessenberg(A, b, x0=x0)
        this.A = A
        x0 = x0 === nothing ? b : x0
        this.x = x0
        this.b = b
        this.residual = b - A(x0)
        return this
    end

    function ConjugateGradient(A::Matrix, b::Union{Matrix, Vector})
        this = ConjugateGradient((x)->A*x, b, x0=b)
        return this
    end

end


# Operator Override
function (this::ConjugateGradient)()
    d, _ = this.ih()
    Ad = this.A(d)
    α = dot(this.residual, d)/dot(d, Ad)
    this.x += α*d
    this.residual = this.b - this.A(this.x)
    # 2 updates to keep the conjugations
    return this.x, sqrt(dot(this.residual, this.residual))
end

