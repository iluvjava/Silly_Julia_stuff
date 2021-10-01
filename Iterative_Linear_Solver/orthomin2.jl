### Implementing the orthomin 2 algorithm for abstract type of linear operators. 
###     * It should converge for all types of linear operators that doesn't have 
###       zero in its field of values. 

import Logging
using LinearAlgebra


function OrthoMin2(A, b, x0; maxitr::Int64=1000, ϵ=1e-3)
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
        if a == 0
            Logging.@warn "Field of values contains zero. "
        end
        x[KK] = x[KK - 1] + a*p[KK - 1]
        r[KK] = r[KK - 1] - a*s[KK - 1]
        p[KK] = s[KK - 1]
        s[KK] = A(s[KK - 1])
        for LL ∈ min(2, KK)
            β = proj(s[KK], s[KK - LL])
            p[KK] = p[KK] - β*p[KK - LL]
            s[KK] = s[KK] - β*s[KK - LL]
        end
        # println("Residual: $(sqrt(abs(dot(r[KK], r[KK]))))")
        if sqrt(abs(dot(r[KK], r[KK]))) < ϵ
            return x[KK]
        end
        
    end
    Logging.@warn "Convergence failed"
    return x[maxitr]
end


