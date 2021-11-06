using Test
using Logging
include("iterative_conjugate_original.jl")

@testset "Testing the original Conjugate Gradient" begin

    function Test1(N=4)
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        IcgOri = IterativeCGOriginal(A, b)   
        
        for _ = 1:N + 1
            display(norm(IcgOri()))    
        end
        @info "Asserting Error is small."
        @assert norm(b - A*IcgOri.x) < 1e-10
        return true
    end

    @test Test1()
    
end