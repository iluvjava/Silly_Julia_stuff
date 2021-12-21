include("CG_Convergence.jl")

# ------------------------------------------------------------------------------
# Interactive Parts: 
# Plotting out the distribution of the Eigenvals. 

n = 30
λ1 = 0.001
fig = plot(title="Eigenvalue Distribution")
Eigenvalues = nothing
for ρ = [0.1, 0.5, 0.8, 0.9, 1]
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
n = 50
fig3 = plot(title="Resnorm2, n=$(n)")
fig4 = plot(title="Error Energy Norm")

for ρ = [0.1, 0.5, 0.8, 0.9, 1]
    A = GetNastyPSDMatrix(ρ, n)  
    b = rand(size(A, 2))
    A = convert(Matrix{Rational{BigInt}}, A)
    b = convert(Vector{Rational{BigInt}}, b)
    cg, Guesses, ResNorm = RunCGTillEnd(
        A, b, 
        maxitr=50, 
        epsilon=1e-10, 
        cg_implementation=Sproj.IterativeCGOriginal
    )
    plot!(fig3, ResNorm, label="ρ=$(ρ)")

    ResEnergyNorm = EnergyErrorNorm(A, b, Guesses)
    plot!(
        fig4, 
        ResEnergyNorm[ResEnergyNorm .!= 0].|>log10,
        label="ρ=$(ρ)"
    )
    
end
display(fig3)
display(fig4)
savefig(fig3, "plots/Exact_Arithematic_Res2Norm.png")
savefig(fig4, "plots/Exact_Arithematic_ResEnergyNorm.png")


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


# ------------------------------------------------------------------------------
# Now let's bound this with the theory we developed. 


