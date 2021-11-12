include("iterative_conjugate_gradient_via_lanczos.jl")
using Test 
using Logging

@testset "Testing the Iterative CG via Lanczos" begin
    function Test1(N=50)
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        cg = IterativeCGViaLanczos(A, b)
        Error = Inf
        MaxItr = 20*N
        while Error > 1e-8 && MaxItr < 5000
            Error = cg()
            println("2 Norm Error: $(Error)")
            MaxItr += 1
        end
        return true
    end
    @test Test1()
end



