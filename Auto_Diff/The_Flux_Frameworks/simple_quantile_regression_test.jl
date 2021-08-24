import BenchmarkTools as Bench

include("simple_quantile_regression.jl")
featureIn = 100
featureOut = 1
N = 100
A = rand(N, featureIn)
y = A*ones(featureIn) + 0.1*rand(N)

quantile = QuantileRegression(featureIn, featureOut, q=0.1)
ŷ = Predict(quantile, A)
l = Loss(quantile, A, y)

losses = []
η=0.01
Bench.@time for II ∈ 1:1000
    for (p, g) ∈ pairs(∇(quantile, A, y))
        p -= η*g
    end
    push!(losses, Loss(quantile, A, y))
end

