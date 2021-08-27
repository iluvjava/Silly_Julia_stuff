include("prototype.jl")

inExactList = [InexactFloat(a)*2^54 for a ∈ rand(2^10)]
exactList = [Rational{BigInt}(a.x) for a ∈ inExactList]
sum1 = sum(inExactList)
sum2 = sum(exactList)
println("Inexact Number yeild: $(sum1)")
println("Exact Rational is: $(sum2)")
println("difference is: $(abs(sum1.x - float(sum2)))")


f(x) = x^2 + 2x + 1
Δx = 1e-10
∂f(x) = (f(x + Δx) - f(x))/Δx
println(∂f(InexactFloat(2.0)))
