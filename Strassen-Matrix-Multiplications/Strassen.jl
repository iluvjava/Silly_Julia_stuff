
module Strassen

    function NaiveMultiply(A:: Array{Float64, 2}, B:: Array{Float64, 2})
        C = Array{Float64}(undef, size(A, 1), size(B, 2))
        for II = 1: size(A, 1)
            for JJ = 1: size(B, 2)
                for KK = 1: size(B, 1)
                    C[II, JJ] += A[II, KK]*B[KK, JJ]
                end
            end
        end
        return C
    end

    function StrassenMultiply(A:: Array{Float64}, B::Array{Float64}; blockSize = 8)

        if (blockSize ≤ 1)
            throw("Blocksize has to be an integer larger than 1.")
        end
        
        α, β = size(A) 
        _, κ = size(B)
        BlockSize = blockSize
        
        # ----------------------------------------------------------------
        # prepare a array of tensors for intermediate storage
        function GetRecursiveDepth(A, B; depth = 1, M = nothing)
            if M === nothing
                M = Array{Array{Float64}}(undef, 0)
            end
            m = size(A, 1) ÷ 2
            n = size(B, 1) ÷ 2
            k = size(B, 2) ÷ 2
            append!(M, [zeros(m, k, 7)])
            if 2m <= BlockSize || 2n <= BlockSize || 2k <= BlockSize
                return depth, M
            end
            return GetRecursiveDepth(view(A, 1:m, 1:n), view(B, 1:n, 1:k); depth=depth + 1, M=M)
        end
        Depth, Mem = GetRecursiveDepth(A, B)
        C = zeros(α, κ)
        
        # ----------------------------------------------------------------
        # Inner multable, recursive func
        function InnerRecursive!(A::Union{SubArray{Float64}, Array{Float64}},
                                B::Union{SubArray{Float64}, Array{Float64}},
                                C::Union{SubArray{Float64}, Array{Float64}};
                                depth:: Int64 = 1)
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
            M = Mem[depth]  # Still a reference copy 
            InnerRecursive!(A11 + A22, B11 + B22, view(M, 1:m, 1:k, 1); depth=depth + 1)
            InnerRecursive!(A21 + A22, B11,       view(M, 1:m, 1:k, 2); depth=depth + 1)
            InnerRecursive!(A11,       B12 - B22, view(M, 1:m, 1:k, 3); depth=depth + 1)
            InnerRecursive!(A22,       B21 - B11, view(M, 1:m, 1:k, 4); depth=depth + 1)
            InnerRecursive!(A11 + A12, B22,       view(M, 1:m, 1:k, 5); depth=depth + 1)
            InnerRecursive!(A21 - A11, B11 + B12, view(M, 1:m, 1:k, 6); depth=depth + 1)
            InnerRecursive!(A12 - A22, B21 + B22, view(M, 1:m, 1:k, 7); depth=depth + 1)

            # Compute and assigned results back, values are garbage collected. 
            C[1: m, 1: k]           = view(M, 1:m, 1:k, 1) + view(M, 1:m, 1:k, 4) - view(M, 1:m, 1:k, 5) + view(M, 1:m, 1:k, 7)
            C[1: m, k + 1: 2k]      = view(M, 1:m, 1:k, 3) + view(M, 1:m, 1:k, 5)
            C[m + 1: 2m, 1: k]      = view(M, 1:m, 1:k, 2) + view(M, 1:m, 1:k, 4)
            C[m + 1: 2m, k + 1: 2k] = view(M, 1:m, 1:k, 1) - view(M, 1:m, 1:k, 2) + view(M, 1:m, 1:k, 3) + view(M, 1:m, 1:k, 6) 
            
            # Handling the Residual rows when extra rows or columns cannot be divisible by 2.
            Bool1 = size(A, 1) % 2 == 1
            Bool2 = size(A, 2) % 2 == 1
            Bool3 = size(B, 1) % 2 == 1
            Bool4 = size(B, 2) % 2 == 1
            x = Bool1 && Bool2 ? A[end, end] : 0
            y = Bool3 && Bool4 ? B[end, end] : 0
            a = Bool1 ? view(A, 2m + 1, 1: 2n) : zeros(2n)
            b = Bool2 ? view(A, 1: 2m, 2n + 1) : zeros(2m)
            d = Bool3 ? view(B, 2n + 1, 1: 2k) : zeros(2k)
            c = Bool4 ? view(B, 1: 2n, 2k + 1) : zeros(2n)
            C[1: 2m, 1: 2k] += b*d'
            C[[end], 1: 2k] += a'*view(B, 1: 2n, 1: 2k) + x.*d'
            C[1: 2m, end]   += view(A, 1: 2m, 1: 2n)*c + y.*b
            C[end, end]     += a'*c + x*y
        end
        # ----------------------------------------------------------------
        InnerRecursive!(A, B, C)
        return C
    end 

    export StrassenMultiply, NaiveMultiply
end 