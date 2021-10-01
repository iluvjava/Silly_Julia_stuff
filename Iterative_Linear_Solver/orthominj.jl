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
    if this.itr == 0 # initialization
        
    end
end




