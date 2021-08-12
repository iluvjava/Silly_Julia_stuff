using COSMO, JuMP, LinearAlgebra, SCS

n = 200
A = randn(n, n)
b = rand(n, 1)
λ = 0.0001

# let's do a 2 norm optimiztion. 
model = Model(with_optimizer(COSMO.Optimizer))
@variable(model, x[1:n])
setvalue.(x, A\b)
@variable(model, η[1:n])
@constraint(model, -η .<= x)
@constraint(model, x .<= η)
@objective(model, Min, λ*sum(η) + sum((A*x - b).^2)) 

optimize!(model)