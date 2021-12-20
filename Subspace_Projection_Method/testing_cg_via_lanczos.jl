include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods
using Test 
using Logging

@testset "Testing the Iterative CG via Lanczos" begin
    function Test1(N=50)
        @info "This is a very basic test to see if the algorithm terminates 
        with reasonable number of steps, from the residual it returns. "
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        cg = Sproj.IterativeCGViaLanczos(A, b)
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

    function Test2(N=50)
        @info "Testing the vadility of the solutions, terminating using the residual returned by the algorithm. "
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        cg = Sproj.IterativeCGViaLanczos(A, b)
        Error = Inf
        MaxItr = 20*N
        Itr = 0
        while Error > 1e-16 && Itr < MaxItr
            Error = cg()
            println("2 Norm Error: $(Error)")
            Itr += 1
        end
        println("Iterations Underwent: $(Itr)")
        ResNorm = norm(b - A*cg.x)
        println("The recomputed residual norm2 is: $(ResNorm)")
        println("The residual informed by the iterative cgl is: $(cg.r)")
        return ResNorm < 1e-6
    end
    function Test3(N=50)
        @info "Testing how well this thing handles exact arithematic."
        A = convert(Matrix{Rational{BigInt}}, rand(N, N))
        b = convert(Vector{Rational{BigInt}}, rand(N))
        A = A*A'
        cg = Sproj.IterativeCGViaLanczos(A, b)
        for Iteration in 1: N + 1
            println(convert(Float64, cg()))
        end
        FinalResNorm = convert(Float64, norm(b - A*cg.x))
        println("Recomputing the residual Norm: $(FinalResNorm)")
        return FinalResNorm < 1e-16
        
    end

    @test Test1()
    @test Test2()
    @test Test3()
end



