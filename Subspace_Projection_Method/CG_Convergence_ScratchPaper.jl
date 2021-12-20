# Test out some of the behaviors of conjugate Gradient under floating point 
# predicted by professor greenbaum. 
# Where is what we want to investigate: 
#   * The convergence of the Conjugate Gradient Method
#   * What are the Eigenvalues of the Tridiagonal Matrix under Finite Precisions 
#   arithemetic compare to the Eigenvalues of the original Matrix? 

include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods
using LinearAlgebra
using Logging
using Plot

"""
    Bad eigen value distribution, they are distributed from 0.001 to 1 like a 
    geometric series. 
    *  A lot of small eigenvalues are clustered close to zero, few larger ones are 
    far away from the other Eigenvalues. 
"""
function GetNastyPSDMatrix(rho::Number, N=20)
    @assert rho <= 1 && rho >= 0
    A = rand(N, N)
    Q, _ = qr(A)
    EigenValues = zeros(N)
    EigenMin, EigenMax = 0.001, 1
    EigenValues[1], EigenValues[end] = EigenMin, EigenMax
    for IdexI in 2:N-1
        EigenValues[IdexI] = EigenMin + 
        ((IdexI - 1)/(N - 1))*(EigenMax - EigenMin)*rho^(N - IdexI)
    end
    return Q*diagm(EigenValues)*Q'
end

function RunCGTillEnd(A, b, Maxitr=100)
    cg = Sproj.IterativeCGViaLanczos(A, b)
    Xs = Vector{typeof(b)}()
    push!(Xs, cg.x)
    Counter = 0
    while cg() < 1e-10 && Counter < Maxitr
        push!(Xs, cg.x)
        Counter += 1
    end
    return cg, Xs
end

function Run()
    TestMatrices = Vector{Matrix{Float64}}()
    
    for rho in 0.1:0.1:0.9
        TestMatrix = GetNastyPSDMatrix(rho)
        b = rand(size(TestMatrices, 1))
        
        push!(TestMatrices, TestMatrix)
    end

    function ProduceConvergenceRate()

    end
    function ProduceEigenValuesPlots()

    end

end


