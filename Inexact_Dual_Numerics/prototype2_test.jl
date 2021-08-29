include("prototype2.jl")

for II in 1:1000
    inExactList = [InexactFloat(a)*2^56 for a ∈ rand(2^8)]
    exactList = [Rational{BigInt}(a.x) for a ∈ inExactList]
    sum1 = sum(inExactList)
    sum2 = sum(exactList)
    println("Inexact Number yeild: $(sum1)")
    println("Exact Rational is: $(sum2)")
    println("difference is: $(abs(sum1.x - float(sum2)))")
    println("Is it in the error bound: $(float(sum2) ∈ sum1)")
    @assert float(sum2) ∈ sum1 "Ok one instance broke. "
end


f(x) = x^2 + 2x + 1
Δx = 1e-15
δf(x) = (f(x - Δx) - f(x))/Δx
println(δf(InexactFloat(4.0)))
Δx = 1e-10
δf(x) = (f(x - Δx) - f(x))/Δx
println(δf(InexactFloat(2.0)))
