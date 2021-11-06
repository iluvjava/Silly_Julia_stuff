include("iterative_conjugate_gradient_via_arnoldi.jl")

using Test

@testset "Tests for Conjugate Gradient Method via Arnoldi" begin
    
    function Test1(N=5)
        @info "Test Run Conjugate Gradient"
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        cg = IterativeCGViaArnoldi(A, b)
        Q = similar(A)
        cg()
        return true
    end

    @test Test1()

end