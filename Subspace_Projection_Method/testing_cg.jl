include("conjugate_gradient.jl")

using Test

@testset "Tests for Conjugate Gradient Method" begin
    
    function Test1(N=5)
        @info "Test Run Conjugate Gradient"
        A = rand(N, N)
        A = A*A'
        b = rand(N)
        cg = ConjugateGradient(A, b)
        Q = similar(A)
        for IdxI in 1:100
            x, res = cg()
            println(norm(b - A*x))
        end 
        return true
    end

    @test Test1()

end