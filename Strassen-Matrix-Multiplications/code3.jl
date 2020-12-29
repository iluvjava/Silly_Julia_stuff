import Base.Threads.@spawn
function StrassenMultiply(A:: Array{Float64}, B::Array{Float64})
    α, β = size(A) 
    _, κ = size(B)
    BlockSize = 128
    # ----------------------------------------------------------------
    # For intermediate storage. 
    # M    = zeros(α, κ, 7)  # Not thread safe. 
    
    C = zeros(α, κ)
    # ----------------------------------------------------------------
    # Inner func
    function InnerRecursive!(A::Union{SubArray{Float64}, Array{Float64}},
                             B::Union{SubArray{Float64}, Array{Float64}},
                             C::Union{SubArray{Float64}, Array{Float64}})
        if size(C, 1) != size(A, 1) || size(C, 2) != size(B, 2) || size(A, 2) != size(B, 1)
            println("Size of A: ", size(A))
            println("Size of B: ", size(B))
            println("Sizs of C: ", size(C))
            throw("Dimension mismatched.")
        end

        m = size(A, 1) ÷ 2
        n = size(B, 1) ÷ 2
        k = size(B, 2) ÷ 2
        if 2m <= BlockSize || 2n <= BlockSize || 2k <= BlockSize
            # Base case here. 
            C[:, :] = A*B
            return 
        end

        # Partition A, reference copy
        A11 = view(A, 1: m, 1: n)
        A12 = view(A, 1: m, n + 1: 2n)
        A21 = view(A, m + 1: 2m, 1: n)
        A22 = view(A, m + 1: 2m, n + 1: 2n)
        
        # Partition B, reference copy
        B11 = view(B, 1: n, 1: k)
        B12 = view(B, 1: n, k + 1: 2k)
        B21 = view(B, n + 1: 2n, 1: k)
        B22 = view(B, n + 1: 2n, k + 1: 2k)
        
        # Freaking Recursion 
        M = zeros(m, k, 7)  # Garbage collected
        sub1 = @spawn InnerRecursive!(A11 + A22, B11 + B22, view(M, 1:m, 1:k, 1))
        sub2 = @spawn InnerRecursive!(A21 + A22, B11,       view(M, 1:m, 1:k, 2))
        sub3 = @spawn InnerRecursive!(A11,       B12 - B22, view(M, 1:m, 1:k, 3))
        sub4 = @spawn InnerRecursive!(A22,       B21 - B11, view(M, 1:m, 1:k, 4))
        sub5 = @spawn InnerRecursive!(A11 + A12, B22,       view(M, 1:m, 1:k, 5))
        sub6 = @spawn InnerRecursive!(A21 - A11, B11 + B12, view(M, 1:m, 1:k, 6))
        sub7 = InnerRecursive!(A12 - A22, B21 + B22, view(M, 1:m, 1:k, 7))
        fetch(sub1)
        fetch(sub2) 
        fetch(sub3) 
        fetch(sub4)
        fetch(sub5)
        fetch(sub6)
        fetch(sub7)

        # Compute and assigned results back, values are garbage collected. 
        C[1: m, 1: k]           = view(M, 1:m, 1:k, 1) + view(M, 1:m, 1:k, 4) - view(M, 1:m, 1:k, 5) + view(M, 1:m, 1:k, 7)
        C[1: m, k + 1: 2k]      = view(M, 1:m, 1:k, 3) + view(M, 1:m, 1:k, 5)
        C[m + 1: 2m, 1: k]      = view(M, 1:m, 1:k, 2) + view(M, 1:m, 1:k, 4)
        C[m + 1: 2m, k + 1: 2k] = view(M, 1:m, 1:k, 1) - view(M, 1:m, 1:k, 2) + view(M, 1:m, 1:k, 3) + view(M, 1:m, 1:k, 6) 
        M = Nothing
        # Handling the Residual rows


    end

    function GetRecursiveDepth(A, B, depth)
        m = size(A, 1) ÷ 2
        n = size(B, 1) ÷ 2
        k = size(B, 2) ÷ 2
        if m <= BlockSize || n <= BlockSize || k <= BlockSize
            return depth
        end
        return GetRecursiveDepth(view(A, 1:m, 1:n), view(B, 1:m, 1:n), depth + 1)
    end
    # ----------------------------------------------------------------
    println("Recur depth is: ", GetRecursiveDepth(A, B, 1))
    InnerRecursive!(A, B, C)
    return C
end


M1 = rand(2^10, 2^10); M2 = rand(2^10, 2^10)
@time C = StrassenMultiply(M1, M2)
BooleanMatrix = map((x)-> abs(x) ≤ 1e-9, C - M1*M2)
println("")
display(sum(BooleanMatrix))
@time M1*M2