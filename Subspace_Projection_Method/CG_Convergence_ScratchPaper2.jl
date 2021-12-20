# Load up the things first. 
include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods
using LinearAlgebra
using Logging
using Plots


"""
    Bad eigen value distribution, they are distributed from 0.001 to 1 like a 
    geometric series. 
     
    * A lot of small eigenvalues are clustered close to zero, few larger ones are 
    far away from the other Eigenvalues. 

    Parameters:
        rho: Number: 
            An number between 0 and 1 that paramaterize the distribution of the eigenvalues for the PSD matrix. 
        N=20: 
            The size of the PSD matrix.
    
    returns: 
        Q^TΛQ: Where Q is a random unitary matrix. 
"""
function GetNastyPSDMatrix(rho::Number, N=20)
    @assert rho <= 1 && rho >= 0
    A = rand(N, N)
    Q, _ = qr(A)
    EigenValues = zeros(N)
    EigenMin, EigenMax = 0.001, 1    # Min Max Eigenvalues. 
    EigenValues[1], EigenValues[end] = EigenMin, EigenMax
    for IdexI in 2:N-1
        EigenValues[IdexI] = EigenMin + 
            ((IdexI - 1)/(N - 1))*(EigenMax - EigenMin)*rho^(N - IdexI)  # formulas
    end
    return Q*diagm(EigenValues)*Q'
end

"""
    Perform the Conjugate Gradient method using the IterativeCGViaLanczos method. 
    
    Parameters: 
        A: The matrix.
        b: The vector on the RHS. 
    returns

"""
function RunCGTillEnd(A, b, Maxitr=100, epsilon=1e-8)
    cg = Sproj.IterativeCGViaLanczos(copy(A), copy(b))
    Xs = Vector{typeof(b)}()    # list of solutions
    Rs = Vector{Float64}()      # list of 2 norm of residuals, recomputed.  
    push!(Xs, cg.x)             # intial guess added. 
    Counter = 0
    ResNorm = Inf
    
    while ResNorm > epsilon && Counter < Maxitr
        _ = cg()
        push!(Xs, cg.x)
        ResNorm = norm(b - A*cg.x)
        push!(Rs, ResNorm)
        Counter += 1
    end
    return cg, Xs, Rs
end

"""

"""
function EnergyErrorNorm()

end


A = GetNastyPSDMatrix(0.9, 50)  # A random matrix
b = rand(size(A, 2))
cg, Guesses, ResNorm = RunCGTillEnd(A, b)
# Residual Norm Progression: 
fig = plot(title="Resnorm2")
plot!(fig, ResNorm)



