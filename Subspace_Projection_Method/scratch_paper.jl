include("iterative_hessenberg.jl")

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

ih = core()