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
    
    prog = ProgressThresh(epsilon, "Minimizing:")

    while ResNorm > epsilon && Counter < maxitr
        _ = cg()
        push!(Xs, cg.x)
        ResNorm = norm(b - A*cg.x)
        push!(Rs, convert(Float64, ResNorm))
        update!(prog, ResNorm)
        Counter += 1
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


# ------------------------------------------------------------------------------
# Interactive Parts: 
# Plotting out the distribution of the Eigenvals. 

n = 30
λ1 = 0.001
fig = plot(title="Eigenvalue Distribution")
Eigenvalues = nothing
for ρ = [1, 0.9, 0.4]
    EigenValues = Vector{Float64}()
    for i in 1: n - 1 
        push!(EigenValues, λ1 + ((i - 1)/(n - 1))*ρ^(n - i))
    end
    plot!(fig, EigenValues, label="ρ = $(ρ)")
    
end
display(fig)
savefig(fig, "plots/Eigenvalues_Distribution.png")

# ------------------------------------------------------------------------------
# Floating Points Arithematic

n = 50
fig1 = plot(title="Resnorm2, n=$(n)")
fig2 = plot(title="Error Energy Norm")

for ρ = [0.1, 0.5, 0.8, 0.9, 1]
    A = GetNastyPSDMatrix(ρ, n)                                                  # A random matrix
    b = rand(size(A, 2))
    cg, Guesses, ResNorm = RunCGTillEnd(A, b, maxitr=2000, epsilon=1e-10)        # Reassign gloabal.
    plot!(fig1, ResNorm, label="ρ=$(ρ)")
    plot!(fig2, EnergyErrorNorm(A, b, Guesses), yaxis=:log10, label="ρ=$(ρ)")
    
end
display(fig1)
display(fig2)
savefig(fig, "plots/Floats_Convergence_Resnorm2.png")
savefig(fig, "plots/Floats_Convergence_ResEnergyNorm.png")

# ------------------------------------------------------------------------------
# Exact Arithematic for the Lanczos

fig3 = plot(title="Resnorm2, n=$(n)")
fig4 = plot(title="Error Energy Norm")

for ρ = [0.1, 0.5, 0.8, 0.9, 1]
    A = GetNastyPSDMatrix(ρ, n)  
    b = rand(size(A, 2))
    A = convert(Matrix{Rational{BigInt}}, A)
    b = convert(Vector{Rational{BigInt}}, b)
    cg, Guesses, ResNorm = RunCGTillEnd(
        A, b, 
        maxitr=2000, 
        epsilon=1e-10, 
        cg_implementation=Sproj.IterativeCGOriginal
    )
    plot!(fig3, ResNorm, label="ρ=$(ρ)")
    plot!(
        fig4, 
        EnergyErrorNorm(A, b, Guesses) .|> log10.|> InfNan2Zero, 
        label="ρ=$(ρ)"
    )
    
end
display(fig3)
display(fig4)
savefig(fig3, "code/Exact_Arithematic_Res2Norm.png")
savefig(fig4, "code/Exact_Arithematic_ResEnergyNorm.png")


# ------------------------------------------------------------------------------
# Get the T matrix for some set value of rho.

ρ = 0.8
n = 30     # !!! n is changed here. 
T = nothing
TrueEigenVals = nothing
cg = nothing
let 
    A = GetNastyPSDMatrix(ρ, n)                                                  # A random matrix
    b = rand(size(A, 2))
    A = convert(Matrix{Float16}, A)
    b = convert(Vector{Float16}, b)
    LocalCG, Guesses, ResNorm = RunCGTillEnd(A, b, maxitr=200, epsilon=1e-14)
    global T = Sproj.GetTMatrix(LocalCG.il)
    global TrueEigenVals, _ = eigen(A)
    global cg = LocalCG
    return 
end

fig5 = scatter(zeros(length(TrueEigenVals)), TrueEigenVals.|>log10, size=(1000,2000), dpi=250)

v = nothing
for Indx in 3:size(T, 1)
    Jndx = Indx - 2
    v, _ = eigen(convert(Matrix{Float16}, T[1:Indx, 1:Indx]))
    scatter!(fig5, zeros(length(v)) .+ Jndx, v.|>log10, legend=false, markershape=:cross)

end


display(fig5)
savefig(fig5, "plots/Lanczos_EigenValues_Convergence.png")
@info "EigenValues from lanczos"
display(v)
display(TrueEigenVals)

# Save the intermediate results for analysis. 
CSV.write("data/T.csv", CSV.Tables.table(T))
CSV.write("data/TrueEigenvalues.csv", CSV.Tables.table(TrueEigenVals))