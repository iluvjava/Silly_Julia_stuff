
using LinearAlgebra

function ArnoldiIterate!(A::AbstractArray{T}, Q::Matrix{T}, H::Matrix{T}) where {T<: ComplexF64}
"""
    prototype first. 
#Arguments
- `A`: A complex matrix

Returns Q, H, where Q is the unitary matrix and H is in 
upper Hessenberg Form. 
"""
    @assert ndims(A) == 2 "A needs to be a 2d matrix"
    m, n = size(A)
    @assert m == n "A needs to be a square matrix"
    @assert size(A) == size(Q) == size(H) "A, Q, H should have the same size"
    q1 = rand(n)  + rand(n)*im 
    q1 ./= norm(q1)
    Q[:, 1] .= q1
    q = similar(q1)
    for k in 2:n + 1
        q .= A*Q[:, k - 1]
        # can be optimized?
        for l in 1:k - 1
            H[l, [k - 1]] .= Q[:, l]'*q
            q -= H[l, k - 1].* Q[:, l]
        end
        # H[1:k - 1, k - 1] = Q[:, 1:k - 1]'*q
        # q -= Q[:, 1: k - 1]*H[1:k - 1, k - 1]
        if k <= n
            H[k, k - 1] = norm(q)
            Q[:, k] = q./norm(q)
        end
    end
    return Q, H
end

