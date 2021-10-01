### Implementing the orthomin 2 algorithm for abstract type of linear operators. 
###     * It should converge for all types of linear operators that doesn't have 
###       zero in its field of values. 

import Logging
using LinearAlgebra
using Test

function OrthoMin2(A, b, x0; maxitr::Int64=1000, ϵ=1e-10, verbose::Bool=false)
    """
        implementation
    
    """
    function proj(x, y)
        return dot(x, y)/dot(y, y)
    end
    r, p, s, x = Dict(), Dict(), Dict(), Dict() # all variables
    r[0] = b - A(x0)
    p[0] = r[0]
    s[0] = A(p[0])
    x[0] = x0
    for KK ∈ 1:maxitr
        a = proj(r[KK - 1], s[KK - 1])
        if a ≈ 0
            Logging.@warn "Field of values contains zero. I refuse to compute"
            throw(ErrorException("0 is in the field of value for the given linear operator"))
        end
        x[KK] = x[KK - 1] + a*p[KK - 1]
        r[KK] = r[KK - 1] - a*s[KK - 1]
        # r[KK] = b - A(x[KK])  # TODO: Investigate this alternative. 
        p[KK] = s[KK - 1]
        s[KK] = A(s[KK - 1])
        for LL ∈ 1:min(2, KK)
            β = proj(s[KK], s[KK - LL])
            p[KK] = p[KK] - β*p[KK - LL]
            s[KK] = s[KK] - β*s[KK - LL]
        end
        NormR = r[KK]⋅r[KK]|>abs|>sqrt
        NormS = s[KK]⋅s[KK]|>abs|>sqrt
        if verbose
            print("Residual: $(convert(Float64, NormR)); ")
            println("s-norm: $(convert(Float64,NormS));")
        end
        if isnan(NormS) || isinf(NormS)
            println("Pause!")
            Logging.@warn("The norm of s is inf or nan and zero is not in the field of value of x")
            throw(ErrorException("Stuff Blowed up."))
        end
        if NormS < ϵ 
            return x[KK]
        end
    end
    Logging.@warn "Convergence failed; maxitr exceeded"
    return x[maxitr]
end

function TestIt1(N)
    A = diagm(rand(N)) .+ 1
    println("Test1, positive diagonal matrix with maxitr $N, square matrix of size $N")
    OrthoMin2((x)-> A*x, ones(N), rand(N); verbose=true, maxitr=2*N)
    println("Finished.")
    println()
    return true
end


function TestIt2(N)
    println("Test2, positive diagonam square matrix with size $N with rational arithmetic")
    A = diagm(rand(N)) .+ 1
    A = convert(Matrix{Rational{BigInt}}, A)
    b = convert(Vector{Rational{BigInt}}, rand(N))
    res = OrthoMin2((x)-> A*x, b, b; verbose=true, maxitr=N)
    display(convert(Vector{Float64}, res))
    return true
end


function TestIt3(N)
    A = rand(N, N)
    A .+= (A|>x->reshape(x, length(x))|>x->sum(x, dims=1)|>diagm)
    println("Test3, random positive with maxitr $N, square matrix of size $N")
    OrthoMin2((x)-> A*x, ones(N), rand(N); verbose=true, maxitr=1000)
    println("Finished.")
    println()
end

@testset "Basic" begin
    @test TestIt1(10)  # works ok... 
    @test TestIt2(10)
    @test TestIt3(10)    
end
