### Orthogonalizing all previous direction using modified gram-schidmtz. 
### And then use them to reduce the residual vector under the 2-norm. 

import Logging
using LinearAlgebra

mutable struct OrthoMin
    A::Function # linear operator
    b
    x0
    function OrthoMin(A::Function, b, x0)
        new(A, b, x0)
    end
end






