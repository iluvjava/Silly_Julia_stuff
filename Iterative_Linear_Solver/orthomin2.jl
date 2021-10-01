### Implementing the orthomin 2 algorithm for abstract type of linear operators. 
###     * It should converge for all types of linear operators that doesn't have 
###       zero in its field of values. 

import Logging
using LinearAlgebra


function OrthoMin2(A, b, x0; maxitr::Int64=1000, ϵ=1e-3, verbose::Bool=false)
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
        # r[KK] = r[KK - 1] - a*s[KK - 1]
        r[KK] = b - A(x[KK])  # TODO: Investigate this alternative. 
        p[KK] = s[KK - 1]
        s[KK] = A(s[KK - 1])
        for LL ∈ min(2, KK)
            β = proj(s[KK], s[KK - LL])
            p[KK] = p[KK] - β*p[KK - LL]
            s[KK] = s[KK] - β*s[KK - LL]
        end
        e = sqrt(abs(dot(r[KK], r[KK]))) 
        if verbose
            print("Residual: $(e); ")
            println("s norm: $(dot(s[KK], s[KK]) |> abs |> sqrt)")
        end
        if isnan(e) || isinf(e)
            println("Pause!")
            Logging.@warn("The residual is nan, or inf, and zero is not in the field of value of x")
            throw(ErrorException("The residual blowed up"))
        end
        if e < ϵ
            return x[KK]
        end
    end
    Logging.@warn "Convergence failed; maxitr exceeded"
    return x[maxitr]
end


function TestIt(N)
    A = diagm(rand(N)) .+ 1 
    OrthoMin2((x)-> A*x, ones(N), rand(N); verbose=true, maxitr=1000)

end

TestIt(10)  # It never works. 