using Test
using UnicodePlots
include("L1_Lasso.jl")

function VanderMonde(x::Vector{T}, deg::Int64) where {T<: Number} ::Matrix{T}
    

end



function Test1()::Bool
    X = rand(4, 4)
    X .-= mean(X, dims=1)
    y = rand(4, 1)
    β, p = LassoRidgeElasticNet(X, y)
    println(p)
    return true
end

@test Test1()


function Test2()::Bool
    # Test with some regression
    f(x) = cos(x)
    x = LinRange(0, 1, 100)

end