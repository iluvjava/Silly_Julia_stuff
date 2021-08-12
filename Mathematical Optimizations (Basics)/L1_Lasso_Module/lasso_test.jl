include("lasso.jl")
include("utils.jl")

using Test
using UnicodePlots

Test1Results = 0

function Test1()
    """
        Test 1, some basic test on a random matrix, copied columns of matrix. 
    """
    deg = 2
    N = 6
    B = rand(N, deg)
    A = fill(0.0, (N, 2deg))
    A[:, 1:deg] = B
    A[:, deg + 1:end] = 2B.^2 + B

    b = A*ones(2deg, 1)
    instance = LassoSCOP(A, b)
    results = LassoPath(instance)
    display(results)
    
    # Let's plot this
    global Test1Results = results

    ylim = (minimum(Test1Results), maximum(Test1Results))
    ThePlot = lineplot(Test1Results[1, :], ylim=ylim, title="Lasso Path")
    for II ∈ 2:size(Test1Results, 1)
        lineplot!(ThePlot, Test1Results[II, :])
    end
    display(ThePlot)
    ThePlot = heatmap(results[end:-1:begin, :])
    display(ThePlot)
    return true

end

@test Test1()


function Test2()
    N = 10
    deg = 10
    xgrid = collect(LinRange(-2π, 2π, N))
    f(x) = sin.(x)
    y = f(xgrid) + 0.1.*rand(N)
    V = VanderMonde(reshape(xgrid, (N, 1)), deg)
    instance = LassoSCOP(V, reshape(y, (N, 1)))
    results = LassoPath(instance)
    display(results)

    ylim = (minimum(results), maximum(results))
    ThePlot = lineplot(results[1, :], ylim=ylim)
    for II ∈ 2:size(results, 1)
        lineplot!(ThePlot, results[II, :])
    end
    display(ThePlot)

    YHat = V*results[:, end]
    ThePlot = lineplot(xgrid, y)
    scatterplot!(ThePlot, xgrid, YHat)
    display(ThePlot)
    display(heatmap(results[end:-1:begin, :]))
    
    return true
end

@test Test2()