using BenchmarkTools
## 
# Speed is asserted through type hints on type ambigous rvalue. 
# @simd macro 
# @inbounds assertion

function innerDot(a::Vector{Float64}, b::Vector{Float64})
    l = min(length(a), length(b))::Int64
    s = 0.0
    @simd for II in 1:l
        @inbounds s += a[II]*b[II]
    end
    return s
end

function innerdotSlow(a, b)
    l = min(length(a), length(b))
    s = 0
    for II in 1:l
        s += a[II]*b[II]
    end
    return s
end
# For compilation
innerDot([1.1, 2,2], [2.2, 2.2])
innerdotSlow([1.1, 2,2], [2.2, 2.2])

GC.gc() 
function timeIt()
    N = 2^10
    reps = 1000
    a, b = rand(N).-0.5, rand(N).-0.5
    println("Random vector Generated ")
    GC.gc()
    s = 0 
    time = @elapsed for II in 1:reps
        s += innerDot(a, b)
    end
    println(string("time (sec):", time))
    println(string("results: ",s, " "))
    Gflops = (N*2*reps*1e-9)/(time)
    println(string("Gflops/s: ", Gflops))
end
timeIt()