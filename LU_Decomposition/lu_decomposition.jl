
function LUDecompose(A::Matrix{T}) where {T<:Number}
    m, n = size(A)
    @assert m == n "Matrix is not squared, can't be LU decomposed. "
    
    L = similar(A)
    U = similar(A)
    LUDecompose!(A, L, U)
    return L, U
end


function LUDecompose!(A::Matrix{T}, L::Matrix{T}, U::Matrix{T}) where {T<:Number}
    m, n = size(A)
    @assert m == n "Matrix is not squared, can't be LU decomposed. "
    @assert size(A) == size(L) "size of A and L should be the same."
    @assert size(A) == size(U) "size of A and U should be the same. "
    U .= A
    for k = 1:n - 1
        L[k, k] = 1
        L[k, k + 1:end] .= 0
        L[k + 1:end, k] = U[k + 1:end, k]/U[k, k]
        U[k + 1:end, k + 1:end] -= L[k + 1:end, k]*U[k, k + 1:end]'
        U[k + 1:end, k] .= 0
    end
    L[end, end] = 1
end

function PLU!(A::Matrix{T}, L::Matrix{T}, U::Matrix{T}, P::Matrix{T}) where {T<:Number}
    m, n = size(A)
    @assert m == n "Matrix is not squared, can't be LU decomposed. "
    @assert size(A) == size(L) "size of A and L should be the same."
    @assert size(A) == size(U) "size of A and U should be the same. "
    @assert size(A) == size(P) "size of A and P should be the same. "
    U .= A
    P .= 0
    for II in 1:size(P, 1)
        P[II, II] = 1
    end
    for k = 1:n - 1
        L[k, k] = 1
        L[k, k + 1:end] .= 0
        m = argmax(abs.(U[k:end, k]))
        P[:, m], P[:, k] = P[:, k], P[:, m]
        U[k, :], U[m, :] = U[m, :], U[k, :]
        L[k + 1:end, k] = U[k + 1:end, k]/U[k, k]
        U[k + 1:end, k + 1:end] -= L[k + 1:end, k]*U[k, k + 1:end]'
        U[k + 1:end, k] .= 0
    end
    L[end, end] = 1
end

function PLU(A::Matrix{<:Number})
    P = similar(A)
    L = similar(A)
    U = similar(A)
    PLU!(A, L, U, P)
    return P, L, U
end