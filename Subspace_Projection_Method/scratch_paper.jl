using Test 
using LinearAlgebra 
using Logging
include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods

function GetNastyPSDMatrix(rho::Number, N=20)
    @assert rho <= 1 && rho >= 0
    A = rand(N, N)
    # Q, _ = qr(A)
    EigenValues = zeros(N)
    EigenMin, EigenMax = 0.001, 1    # Min Max Eigenvalues. 
    EigenValues[1] = EigenMin
    for IdexI in 2:N
        EigenValues[IdexI] = EigenMin + 
            ((IdexI - 1)/(N - 1))*(EigenMax - EigenMin)*rho^(N - IdexI)  # formulas
    end
    return diagm(EigenValues)
end

A = GetNastyPSDMatrix(0.1)


# Nothing wrong with Lanzos. 