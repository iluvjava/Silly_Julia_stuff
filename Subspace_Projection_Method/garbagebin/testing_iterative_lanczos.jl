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
    
    function Test2(n=5)
        @info "Verifying the recurrence of T, Q, and A matrix."
        A = rand(n, n)
        A = A'+ A # symm
        b = rand(n)
        il = IterativeLanczos(A, b)
        betas = il(n - 1)
        qs = GetPrevious3OrthogonalVec(il)
        Q = GetQMatrix(il)
        T = GetTMatrix(il)
        println("Q Matrix: ")
        display(Q)
        println("T matrix:") 
        display(T)
        # Check error
        recurrError = A*Q[:, 1: end - 1] - Q[:, 1: end - 1]*T -
            betas[end]*(Q[:, end]*Matrix(I, size(Q, 2) - 1, size(Q, 2) - 1)[:, end]')
        recurrError = norm(recurrError, Inf)
        inTol1 = recurrError < 1e-10

        # Check Error
        

        if !inTol1 
            @warn "The recurrence of AQ = QT̃ broke with inf norm of $inTol1, plase check."
        end

        # Check Error
        println("Q^TAQ is:")
        QAQ = Q[:, 1: end - 1]'*A*Q[:, 1: end - 1]
        display(QAQ)
        inTol2 = norm(QAQ - T, Inf) < 1e-10
        if !inTol2
            @warn "Factorizations Q'AQ - T has an inf norm of $inTol2, something broke. "
        end

        return inTol2 && inTol2
    end

    @test Test1()
    @test Test2()
end

