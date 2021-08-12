using Convex, LinearAlgebra
using SCS, COSMO

# initialize data
n = 100
y = rand(n)
x = Variable(n)

# first solve
lambda = 100
problem = minimize(sumsquares(y - x) + lambda * sumsquares(x - 10))
println("=====================================================================")
println("This is the First solve ")
println("=====================================================================")
@time solve!(problem, COSMO.Optimizer)
println("=====================================================================")
println("that was the time it took for the first solve: ")
println("=====================================================================")

# now warmstart
# if the solver takes advantage of warmstarts, 
# this run will be faster
lambda = 100

# problem = minimize(sumsquares(y - x) + lambda * sumsquares(x - 10))
println("=====================================================================")
println("This is the Second solve ")
println("=====================================================================")
@time solve!(problem, COSMO.Optimizer, warmstart=true)
println("=====================================================================")
println("that was the time it took for the second solve: ")
println("=====================================================================")