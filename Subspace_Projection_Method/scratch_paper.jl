using Test 
using LinearAlgebra 
using Logging
include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods

# What is going on with the error produced by the cg via lanczos? 

A = rand(10, 10)
A = A'*A
b = rand(10)
ExactCG = Sproj.IterativeCGOriginal(
    convert(Matrix{Rational{BigInt}}, A), 
    convert(Vector{Rational{BigInt}}, b)
)
LanczosCG = Sproj.IterativeCGViaLanczos(
    A, b
)
for Itr in 1:13
    Resnorm1 = ExactCG()
    Resnorm2 = LanczosCG()
    Resnorm3 = norm(b - A*LanczosCG.x)
    println("Resnorm Exact: $(Resnorm1), Resnorm Lanczos: $(Resnorm2), Resnorm Lanczos Recomp: $(Resnorm3)")
end

