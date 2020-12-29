# This is a failed attempt


function StrassenMultiply(A:: Array{Float64}, B::Array{Float64})
    α, β = size(A) 
    _, κ = size(B)
    # ----------------------------------------------------------------
    # For intermediate storage. 
    M    = zeros(α, κ, 7)  # Not thread safe. 
    C    = zeros(α, κ)
    # ----------------------------------------------------------------
    # Inner func
    function InnerRecursive!(A,
                             B,
                             C
    # Store to this sub-matrix, it's actually the M matrix recursively
                            )
        if size(C, 1) != size(A, 1) || size(C, 2) != size(B, 2) || size(A, 2) != size(B, 1)
            println("Size of A: ", size(A))
            println("Size of B: ", size(B))
            println("Sizs of C: ", size(C))
            throw("Dimension mismatched.")
        end

        m = size(A, 1) ÷ 2
        n = size(B, 1) ÷ 2
        k = size(B, 2) ÷ 2
        if m <= 1 || n <= 1 || k <= 1
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
        M[1:m, 1:k, :] .= 0
        InnerRecursive!(A11 + A22, B11 + B22, view(M, 1:m, 1:k, 1))
        InnerRecursive!(A21 + A22, B11,       view(M, 1:m, 1:k, 2))
        InnerRecursive!(A11,       B12 - B22, view(M, 1:m, 1:k, 3))
        InnerRecursive!(A22,       B21 - B11, view(M, 1:m, 1:k, 4))
        InnerRecursive!(A11 + A12, B22,       view(M, 1:m, 1:k, 5))
        InnerRecursive!(A21 - A11, B11 + B12, view(M, 1:m, 1:k, 6))
        InnerRecursive!(A12 - A22, B21 + B22, view(M, 1:m, 1:k, 7))

        # Compute and assigned results back. 
        C11 = view(M, 1:m, 1:k, 1) + view(M, 1:m, 1:k, 4) - view(M, 1:m, 1:k, 5) + view(M, 1:m, 1:k, 7)
        C12 = view(M, 1:m, 1:k, 3) + view(M, 1:m, 1:k, 5)
        C21 = view(M, 1:m, 1:k, 2) + view(M, 1:m, 1:k, 4)
        C22 = view(M, 1:m, 1:k, 1) - view(M, 1:m, 1:k, 2) + view(M, 1:m, 1:k, 3) + view(M, 1:m, 1:k, 6)

        C[1: m, 1: k]           = C11
        C[1: m, k + 1: 2k]      = C12
        C[m + 1: 2m, 1: k]      = C21
        C[m + 1: 2m, k + 1: 2k] = C22  # Hit point
        # Handling the Residual rows
    end
    # ----------------------------------------------------------------
    InnerRecursive!(A, B, C)
    return C
end


M1 = ones(2^3, 2^3); M2 = ones(2^3, 2^3)
C = StrassenMultiply(M1, M2)
BooleanMatrix = map((x)-> abs(x) ≥ 1e-9, C - M1*M2)
show(stdout, "text/plain", BooleanMatrix)
println("")
display(C)