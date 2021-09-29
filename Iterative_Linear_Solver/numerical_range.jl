# Computing the numerical range of a matrix. 
using LinearAlgebra
using UnicodePlots
import Plots as Plots

function GetSupportLineContactPoint(A::Matrix{Complex}, θ)
    c = exp(im*θ) # rotation 
    H = (c*A + conj(transpose(c*A)))/2  # Hermitian 
    Eigs = eigen(H)
    MaxλθIdx = argmax(Eigs.values)
    Maxμθ = Eigs.vectors[:, MaxλθIdx]
    Results =  dot(Maxμθ, c*A, Maxμθ)
    # println(" θ: $(θ)\n eigs: $(Eigs.values)\n support point, x_θ^HAx_θ: $(Results)")
    return Results
end

function main()
    λ = 0.5
    A = [λ 1; 0 λ]
    CoordsX = Vector{Float64}()
    CoordsY = Vector{Float64}()
    for θ ∈ range(0, 2π, length=900)
        Contact = GetSupportLineContactPoint(convert(Matrix{Complex}, A), θ)
        Contact *= exp(-θ*im)
        # println("θ: $(θ); Contact: $(Contact), absval: $(abs(Contact))")
        push!(CoordsX, real(Contact))
        push!(CoordsY, imag(Contact))
    end
    scatterplot(CoordsX, CoordsY)
    # Plots.scatter(CoordsX, CoordsY)
end

main()
