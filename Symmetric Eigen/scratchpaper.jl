using LinearAlgebra

function PureQRAuto(m::AbstractArray{<:Number}; Abstol:: Number=1e-6)
    if m' != m
        throw("Matrix is not hermitial, don't use PureQRAuto.")
    end 
    MaxItr = 1e4; Counter = 0;
    A = copy(m)
    Q = I
    while norm(A - Diagonal(A), Inf) > Abstol && Counter < MaxItr
        Decomp = qr(A)
        Q     *= Q*Decomp.Q
        A      = Decomp.R*Decomp.Q
        Counter += 1
    end
    print("Error infinity norm: ", norm(A - Diagonal(A), Inf))
    return diag(A), Q
end

function RandTridiognalHermitian(n:: Int64)
    Result = zeros(n, n)
    for II = 1: n
        Result[II, II] = rand()
    end
    RandomArray = rand(1, n - 1)
    for JJ = 1: n - 1
        Result[JJ, JJ + 1] = RandomArray[JJ]
        Result[JJ + 1, JJ] = RandomArray[JJ]
    end
    return Result
end


TridRandMatrix = RandTridiognalHermitian(100)
display(TridRandMatrix)
A, Q = PureQRAuto(TridRandMatrix)
display(A)
