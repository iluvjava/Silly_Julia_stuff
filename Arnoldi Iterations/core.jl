
using LinearAlgebra

function ArnoldiIterate(A)
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
    q1 = rand(n, 1)  + rand(n, 1)*im
    q1 = q1./norm(q1)
    H = similar(A, ComplexF64)
    Q = similar(A, ComplexF64)
    Q[:, 1] = q1
    for k in 2:n + 1
        q = A*Q[:, k - 1]
        for l in 1:k - 1
            H[l, k - 1] = Q[:, l]'*q
            q -= H[l, k - 1].* Q[:, l]
        end
        if k <= n
            H[k, k - 1] = norm(q)
            Q[:, k] = q./norm(q)
        end
    end
    return Q, H
end

