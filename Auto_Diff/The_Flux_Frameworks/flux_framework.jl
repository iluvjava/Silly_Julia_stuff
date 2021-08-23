using Flux
if !isdefined(Base, :A)
    const A = rand(100, 10)
    const y = rand(100)
end

w, b = rand(10), rand(1)

ŷ(w, b) = A*w .+ b
l(y, w, b) = 0.5*sum((ŷ(w, b) - y).^2)

∇l = gradient(()->l(y, w, b), params(w, b))

