include("iterative_lanczos.jl")

using Test
using Logging

@testset "Testing the Iterative Lanczos LDL Algorithm; Basic" begin
    function Test1(N=5)
        A = rand(N, N)
        b = rand(N)
        A = A*A'
        il = IterativeLanczos(A, b)
        il()
        Q = GetQMatrix(il)
        L = GetLMatrix(il)
        T = GetTMatrix(il)
        D = GetDMatrix(il)
        print("QLTD Bass cases: :")
        display(Q)
        display(L)
        display(T)
        display(D)
        for _ in 1:N-1
            R2norm = il()
            Q = GetQMatrix(il)
            L = GetLMatrix(il)
            T = GetTMatrix(il)
            D = GetDMatrix(il)
            println("QLTD mattrices:")
            display(Q)
            display(L)
            display(T)
            display(D)
            Error = norm(reshape(Q'*A*Q - T, :), Inf)
            @assert Error <= 1e-9 "Q'AQ reconstruct error too big: $(Error)"
            Error = norm(reshape(L*D*L' - T, :), Inf)
            @assert Error <= 1e-9 "LDL reconstruct error too big: $(Error)"
            Error = norm(reshape(il.Linv - inv(L)[:, 1], :), Inf)
            @assert Error <= 1e-9 "First column of L^{-1} error too big: $(Error)"
            println("2norm of residual at current step is: $(R2norm)")
        end
        
        return true
    end
    function Test2(N=5)
        A = rand(N, N) + im*rand(N, N)
        b = rand(N) + im*rand(N)
        A = A*A'
        il = IterativeLanczos(A, b)
        il()
        Q = GetQMatrix(il)
        L = GetLMatrix(il)
        T = GetTMatrix(il)
        D = GetDMatrix(il)
        print("QLTD Bass cases: :")
        display(Q)
        display(L)
        display(T)
        display(D)
        for _ in 1:N-1
            R2norm = il()
            Q = GetQMatrix(il)
            L = GetLMatrix(il)
            T = GetTMatrix(il)
            D = GetDMatrix(il)
            println("QLTD mattrices:")
            display(Q)
            display(L)
            display(T)
            display(D)
            Error = norm(reshape(Q'*Q - I, :), Inf)
            @assert Error <= 1e-9 "The orthogonalization error is too big: $(Error)"
            Error = norm(reshape(Q'*A*Q - T, :), Inf)
            @assert Error <= 1e-9 "Q'AQ reconstruct error too big: $(Error)"
            Error = norm(reshape(L*D*L' - T, :), Inf)
            @assert Error <= 1e-9 "LDL reconstruct error too big: $(Error)"
            Error = norm(reshape(il.Linv - inv(L)[:, 1], :), Inf)
            @assert Error <= 1e-9 "First column of L^{-1} error too big: $(Error)"
            println("2norm of residual at current step is: $(R2norm)")
        end
        
        return true
    end
    @info "Testing Real Hermitian"
    @test Test1()
    @info "Testing Complex Hermitian"
    @test Test2()
end