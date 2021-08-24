using Flux, CUDA
import BenchmarkTools as Bench
using LinearAlgebra
m, n = 100,1000
A = rand(m, n) |> gpu
y = rand(m) |> gpu


function start()
    
    w, b = rand(n) |> gpu, rand(1) |> gpu
    function Loss(w, b)
        return sum((A*w .+ b - y).^2)
    end
    losses = []
    η = 1/(2*opnorm(A'*A))
    v = similar(w)
    v .= 0
    @time for II ∈ 1:10000
        loss, ∇ = Flux.Zygote.withgradient(()->Loss(w, b), params(w, b))
        push!(losses, loss)
        v .= 0.9*v - η*∇[w]
        w .+= - η*∇[w]
        b -= η*∇[b]
    end
    display(convert(Vector{Float64}, losses))
    return w, b
end 

w, b = start()
