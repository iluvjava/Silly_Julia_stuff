include("simple_quantile_regression.jl")
featureIn = 10
featureOut = 1
N = 100
A = rand(N, featureIn)
y = A*ones(featureIn) + 0.1*rand(N)

quantile = QuantileRegression(featureIn, featureOut)
ŷ = Predict(quantile, A)
l = Loss(quantile, A, y)

losses = []
η=0.01
for II ∈ 1:100
    ∇l = ∇(quantile, A, y)
    ∂W, ∂b = ∇l
    quantile.W .-= (η)*∂W
    quantile.b .-= (η)*∂b
    push!(losses, Loss(quantile, A, y))
end

η=1e-4
for II ∈ 1:100
    ∇l = ∇(quantile, A, y)
    ∂W, ∂b = ∇l
    quantile.W .-= (η)*∂W
    quantile.b .-= (η)*∂b
    push!(losses, Loss(quantile, A, y))
end