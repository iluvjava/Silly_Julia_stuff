using Test 
using LinearAlgebra 
using Logging
include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods

@testset "Testing the original Conjugate Gradient" begin

    function Test1(N=4)
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        IcgOri = Sproj.IterativeCGOriginal(A, b)
        for _ = 1:N + 1
            display(norm(IcgOri()))    
        end
        @info "Asserting Error is small."
        @assert norm(b - A*IcgOri.x) < 1e-10
        return true
    end

    function Test2(N=10)
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        A = convert(Matrix{Rational{BigInt}}, A)
        b = convert(Vector{Rational{BigInt}}, b)
        IcgOri = Sproj.IterativeCGOriginal(A, b)
        for _ = 1:N + 1
            display(norm(IcgOri()))    
        end
        @info "Asserting Error is small."
        @assert norm(b - A*IcgOri.x) < 1e-10
        return true
    end

    @info "Basic Tests"
    @test Test1() 
    @info "Testing CG original under exact arithematic"
    @test Test2()

end


