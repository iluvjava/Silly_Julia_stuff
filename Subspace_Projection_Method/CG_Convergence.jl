# Load up the things first. 
include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods
using LinearAlgebra
using Logging
using Plots
using Arpack
# gr(); # Specified Backend for Plots
using ProgressMeter
using CSV 


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

"""
    Perform the Conjugate Gradient method using the IterativeCGViaLanczos method. 
    
    Parameters: 
        A: The matrix.
        b: The vector on the RHS. 
    returns

"""
function RunCGTillEnd(A, b; maxitr=100, epsilon=1e-8, cg_implementation=nothing)
    if cg_implementation === nothing
        cg_implementation = Sproj.IterativeCGViaLanczos
    end
    cg = cg_implementation(copy(A), copy(b))
    Xs = Vector{typeof(b)}()    # list of solutions
    Rs = Vector{Float64}()      # list of 2 norm of residuals, recomputed.  
    push!(Xs, cg.x)             # intial guess added. 
    Counter = 0
    ResNorm = Inf
    
    prog = ProgressUnknown("Minimizing Residual, Counter:")

    while ResNorm > epsilon && Counter < maxitr
        _ = cg()
        push!(Xs, cg.x)
        ResNorm = convert(Float64, norm(b - A*cg.x))
        push!(Rs, ResNorm)
        Counter += 1
        update!(prog, Counter)
    end
    finish!(prog) 
    return cg, Xs, Rs
end

"""
    Accept a linear system, and a list of guesses from the Lanczos Algorithm, 
    it will return a list energy norm of the 
    error vector for the linear system. 
    
    Parameters: 
        A: The matrix
        b: The vector on the RHS
        Xs: All the guesses vector from the CGLanczos Algorithm. 
"""
function EnergyErrorNorm(A, b, Xs)
    XStar = A\b
    return Xs .|> (x)-> dot(x - XStar, A, x - XStar)
end

function InfNan2Zero(n)
    if isnan(n) || n == Inf || n == -Inf
        return n
    end
    return n
end
