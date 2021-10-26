include("iterative_hessenberg.jl")

function core(N=3)
    A = rand(N, N)
    A = A*A'
    b = rand(N)
    display(A)
    
    ih = IterativeHessenberg(A, b, max_k=3)
    q1, _ = ih()
    q2, _ = ih()
    q3, _ = ih()
    q4, _ = ih()
    display(ih.H)
    return ih    
end

ih = core()