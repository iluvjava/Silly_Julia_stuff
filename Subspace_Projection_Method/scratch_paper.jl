include("iterative_hessenberg.jl")
import Logging
using SuiteSparse


function core(N=3)
    A = rand(N, N)
    A = A*A'
    b = rand(N)
    x = b
    r = b - A*x
    display(A)
    
    ih = IterativeHessenberg(A, b, max_k=3)
    for _ in 1:4
        q, h = ih()
        α = dot(r, q)/dot(q, A*q)
        x += α*q
        println("α: $(α)")
        r -= α*A*q
        println("new residual norm: $(norm(r)), expect to be $(h)")
        println("new residual dot conjugate dir: $(dot(r, q))")
    end
    
    return
end

function LanczaosLDL()
    function TridiagToSparse(A)
        m, _ = size(A)
        return SymTridiagonal(
            A[[CartesianIndex(I, I) for I in 1:m]],
            A[[CartesianIndex(I, I + 1) for I in 1:m-1]])
    end
    N = 5
    A = rand(N, N)
    A = A*A'
    b = rand(N)
    display(A)
    ih = IterativeHessenberg(A, b, max_k=2)
    ih()
    @info "L Matrices and L Inverses of the Tridiag Matrix"
    for _ in 1:4
        ih()
        T = GetHessenberMatrix(ih)
        Factorization = ldlt(TridiagToSparse(T[1:end-1, :]))
        display(Factorization.L)
        display(inv(Factorization.L))
    end
    
end

function TridiagCheck()
    N = 5
    A = rand(N, N)
    A = A*A'
    b = rand(N)
    display(A)
    ih = IterativeHessenberg(A, b, max_k=2)
    ih()
    @info "Tridiag HESSENBERG" 
    for _ in 1:4
        ih()
        T = GetHessenberMatrix(ih)
        display(T)
    end
    
end

# ih = core()
LanczaosLDL()
TridiagCheck()