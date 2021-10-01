### Orthogonalizing all previous direction using modified gram-schidmtz. 
### And then use them to reduce the residual vector under the 2-norm. 
###   * Store the MGS subspace for the Krylov Subspace 
###   * Support iterators
### Algorithm Features: 
###     Proposed as orthomin(j) methods in Greenbaum's book. It compute krylov 
###     subspace κ and A(κ) on the go and orthogonalizing it on the go. 

import Logging
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

function (this::OrthoMin)(j::Int64=-1)
    A, b, x0 = this.A, this.b, this.x0
    r, x, p, s = this.r, this.x, this.p, this.s
    if this.itr == 0 # initialization
        r[0] = b - A(x0)
        x[0] = x[0]
        p[0] = r[0]
        s[0] = A(p[0])
    end
    for KK ∈ this.itr + 1:this.itr + j
        for LL ∈ 1:KK
            
        end
    end
end




