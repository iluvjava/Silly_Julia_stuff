using Test 
using LinearAlgebra 
using Logging
include("SubspaceProjectionMethods.jl")
Sproj = SubspaceProjectionMethods


# TODO: Investigate why it can't do exact arithematic. 
N = 5
A = rand(N, N)
b = rand(N)
A = A*A'
A = convert(Matrix{Rational{BigInt}}, A)
b = convert(Vector{Rational{BigInt}}, b)
il = Sproj.IterativeLanczos(A, b)

for I in 1:N + 1
    println(il())
end
@info "Matrix L"
display(Sproj.GetLMatrix(il))
@info "Matrix D"
display(Sproj.GetDMatrix(il))
