### Orthogonalizing all previous direction using modified gram-schidmtz. 
### And then use them to reduce the residual vector under the 2-norm. 
###   * Store the MGS subspace for the Krylov Subspace 
###   * Support iterators
### Algorithm Features: 
###     Proposed as orthomin(j) methods in Greenbaum's book. It compute krylov 
###     subspace κ and A(κ) on the go and orthogonalizing it on the go. 

using Logging
using LinearAlgebra

mutable struct OrthoMin
    A::Function # linear operator
    b
    x0
    itr::Int64
    p::Dict
    s::Dict
    x::Dict
    r::Dict
    function OrthoMin(A::Function, b, x0)
        new(A, b, x0, Dict(), 0)
    end
end

function (this::OrthoMin)(j::Int64=1)
    prj(x, y) = x⋅y/y⋅y
    A, b, x0 = this.A, this.b, this.x0
    r, x, p, s = this.r, this.x, this.p, this.s

    if this.itr == 0 # initialization
        r[0] = b - A(x0)
        x[0] = x[0]
        p[0] = r[0]
        s[0] = A(p[0])
    end
    k = this.itr + j
    for KK ∈ this.itr + 1: k
        a = prj(r[KK - 1], s[KK - 1])
        if a ≈ 0 
            @warn "Field of value contains zero. "
            throw(ErrorException("Computational Error"))
        end
        x[KK] = x[KK - 1] + a*p[KK - 1]
        r[KK] = r[KK - 1] - a*s[KK - 1]
        p[KK] = s[KK - 1]
        s[KK] = A(s[KK - 1])
        for LL ∈ 1:KK
            β = proj(s[KK], s[KK - LL])
            p[KK] = p[KK] - β*p[KK - LL]
            s[KK] = s[KK] - β*s[KK - LL]
        end
    end
    this.itr = k
    return r[k]
end

