using Test 
using LinearAlgebra 
using Logging
include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods


# TODO: Investigate this phenomena.``
N = 5
A = rand(N, N)
b = rand(N)
A = A*A'
A = convert(Matrix{Rational{BigInt}}, A)
b = convert(Vector{Rational{BigInt}}, b)
il = Sproj.IterativeLanczos(A, b)
for I in 1:6
    println(il())
end