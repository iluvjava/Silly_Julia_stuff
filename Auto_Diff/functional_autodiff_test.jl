include("functional_autodiff.jl")

# ==============================================================================
# computer
g = Sin^2 + Exp∘(Identity()^2)

# my answer to that derivative
dg(x) = 2sin(x)*cos(x) + 2x*exp(x^2)

println(dg(0.2))
println(g.df(0.2))

# computer
g = Sin∘(2*Identity())

# me
dg(x) = 2*cos(2x)

println(dg(0.5))
println(g.df(0.5))

# comuter
g = Sin∘Exp∘(2*Identity()) 

# me
dg(x) = cos(exp(2x))*(2exp(2x))

println(dg(0.3))
println(∂(g)(0.3))

# ==============================================================================
# f(x) = x, then f(0.3) = 0.3 and ∂f(x) = 1
x = x -> DiffNum(x)
f = x(0.3)
g = 1/(f^2)
println(∂(g))
println(-2/0.3^3)
println(∂(2*x(0.2)^2))
println(4*0.2)
println(∂(x(3)^0.5))
println(0.5*3^(-0.5))
println(∂(exp(x(2)^2)))
println(2*2exp(2^2))

for II in LinRange(-1, 1, 100)
    a = ∂(exp(x(II)^2)*sin(x(II)))
    b = 2II*exp(II^2)*sin(II) + cos(II)*exp(II^2)
    @assert a == b
end

