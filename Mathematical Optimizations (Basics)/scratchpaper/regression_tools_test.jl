using Test
using UnicodePlots

include("regression_tools.jl")


function Test1()::Bool
    X = rand(4, 4)
    X .-= mean(X, dims=1)
    y = rand(4, 1)
    β, p = LassoRidgeElasticNet(X, y)
    println(p)
    println("Coefficnet found is:")
    println(β)
    println("Solving it a second time and test for warm starting: ")
    solve!(p, ()-> SCS.Optimizer(verbose=true), warmstart=true)
    
    return true
end


function Test2()::Bool
    

end




@test Test1()  # Basic solver on the Lasso Ridge Elastic function. 