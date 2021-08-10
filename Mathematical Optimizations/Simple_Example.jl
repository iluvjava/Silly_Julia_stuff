Block1 = begin 
    # cone constraints, convex objective. 
    using Convex, SCS
    m, n = 4, 5
    A = rand(m, n)::Array{Float64}
    b = rand(m, 1)::Array{Float64}
    x = Variable(n)::Variable
    # For floating, relax the constraints, or we will have to delve deep and 
    # tweak the solver. 
    problem = minimize(sumsquares(A*x - b), [x >= 1e-7])::Problem{Float64}
    solve!(problem, () -> SCS.Optimizer(verbose=true))
    problem.status
    display("solutions turns out to be:")
    display(x.value)
end