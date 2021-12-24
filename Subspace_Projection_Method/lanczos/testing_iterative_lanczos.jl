include("iterative_lanczos.jl")
using LinearAlgebra
using Logging
using Test

@testset begin
    function Test1(n=5)
        @info "Testing Lanczos on a $n by $n symmetric matrix real fast. "
        A = rand(n, n)
        A = A' + A # symm
        b = rand(n)
        il = IterativeLanczos(A, b)
        for _ in 1: n - 1
            il()
        end
        println("This is the Q Matrix: ")
        Q = GetQMatrix(il)
        display(Q)
        dotQ = Q'*Q
        inTol = norm(dotQ - I, Inf) < 1e-10
        if !inTol
            @warn "The error deviated from identity turns out to be 
            $(norm(dotQ - I, Inf)), larger than 1e-10, something is wrong. "
        end
        return inTol
    end
    
    @test Test1()

end

