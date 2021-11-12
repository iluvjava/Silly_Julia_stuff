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
        Itr = 0
        while Error > 1e-8 && Itr < MaxItr
            Error = cg()
            println("2 Norm Error: $(Error)")
            Itr += 1
        end
        print("Iterations Underwent: $(Itr)")
        return Itr < MaxItr
    end
    @test Test1()
end



