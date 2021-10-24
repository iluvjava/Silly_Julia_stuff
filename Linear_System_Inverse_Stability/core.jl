using LinearAlgebra
using Logging
# PA = LU in julia. 
# P is the pervec on columns of A 


function MakeProblem(n, ϵ=0.1)
    function HilbertMatrix(n)
        """
            Function for benchmarking algorithms. 
        """
        A = zeros(n, n)
        for IdxI = 1:n, IdxJ = 1:n
            A[IdxI, IdxJ] = 1/(IdxI + IdxJ)
        end
        return A
    end
    A = HilbertMatrix(n) 
    A .+= tril(rand(Float64, size(A)).*ϵ)  # pertubate. 
    x = rand(n, 1)
    Ar = convert(Matrix{Rational{BigInt}}, A)
    xr = convert(Matrix{Rational{BigInt}}, x)
    b = convert(Matrix{Float64}, Ar*xr)
    return A, x, b
end


function BenchTheError(A::T, x::T, b::T) where {T <: AbstractArray}
    Ar = convert(Matrix{Rational{BigInt}}, A)
    xr = convert(Matrix{Rational{BigInt}}, x)
    b̂ = Ar*xr
    b̂ = convert(typeof(b), b̂)
    return norm(b - b̂, 1)/length(b)
end


function JuliaLUDecomp(A::AbstractMatrix)
    L, U, P = lu(A)
    n = size(A)[1]
    P = Matrix{Int}(I, n, n)[P, :]
    return P, L, U
end


function JuliaSolveWithLU(A::T, b::T)where {T<: AbstractArray}
    P, L, U = JuliaLUDecomp(A)
    b = P*b
    b = L\b
    b = U\b
    return b
end

function main()
    print("Testing Hilbert matrix without pertubation")
    

end