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


n = 15
A = convert(Matrix{Float16}, GetNastyPSDMatrix(0.9, n))
b = convert(Vector{Float16}, rand(n))
il = Sproj.IterativeLanczos(A, b)
cg = Sproj.IterativeCGViaLanczos(A, b)
for Idx in 1: 20
    il()
    cg()
end

T1 = Sproj.GetTMatrix(il)
Q1 = Sproj.GetQMatrix(il)
T2 = Sproj.GetTMatrix(cg.il)

# Nothing wrong with Lanzos. 