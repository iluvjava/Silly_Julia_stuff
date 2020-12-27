using LinearAlgebra  # import the library independntly

function SingleShiftedQR!(m::AbstractArray{<: Number}; maxItr::Int64 = 1000, verbose::Bool = true)
    if size(m, 1) != size(m, 2) 
        throw("Matrix not square, cannot run shifted QR algorithm.")
    end
    width = size(m, 1)
    A = m
    μ = A[width, width]
    Counter = 0
    EigenColumn = -1
    QAccumulated = I
    for II = 1: maxItr
        QRDecomp = qr(A - μ.*I)
        Q, R     = QRDecomp.Q, QRDecomp.R
        A        = R*Q + μ.*I
        QAccumulated = QAccumulated*Q
        for JJ = 1: width - 1
            if  abs(A[JJ + 1, JJ]) ≤ 1e-6
                EigenColumn = JJ
                println("EigenColumn = ", EigenColumn)
                @goto OutterLoop
            end
        end 
        Counter += 1
    end
    @label OutterLoop
    println("Number of iterations: ", Counter)
    display(A)
    return A[EigenColumn, EigenColumn], QAccumulated[:, EigenColumn]
end


"""
    Perform a recusive search of the Eigenvalues for the given Hermitian Matrix. 
"""
function RecursiveShiftedQR!(m::AbstractArray{<:Number})
    if size(m, 1) != size(m, 2)
        throw("The shifted QR algorithm is only for square matrices but this is not squared: ")
        display(m)
    end
    if size(m, 1) == 1
        # Base case 
        return
    end
    function SingleShiftedQR!(m::AbstractArray{<: Number}; maxItr::Int64 = 1000)
        A = view(m, :, :)
        width = size(m, 1)
        EigenColumn = -1
        # QAccumulated = I                              # This feature is cancelled in this implementation
        Flag = true
        for II = 1: maxItr
            μ = A[end, end] 
            QRDecomp = qr(A - μ.*I)
            Q, R     = QRDecomp.Q, QRDecomp.R
            A[:, :]  = R*Q + μ.*I
            # QAccumulated = QAccumulated*Q
            for JJ = 1: width - 1
                if  abs(A[JJ + 1, JJ]) ≤ 1e-14
                    EigenColumn = JJ
                    Flag = false
                    @goto OutterLoop
                end
            end 
        end
        @label OutterLoop
        if Flag
            # This really should not happen
            println("Wilskin Shift might have failed to converge.")
        end
        return EigenColumn
    end

    JJ          = SingleShiftedQR!(m)
    AUpper      = view(m, 1: JJ, 1: JJ)
    ALower      = view(m, JJ + 1: size(m, 1),  JJ + 1: size(m, 1))
    RecursiveShiftedQR!(AUpper)  # No return value needed
    RecursiveShiftedQR!(ALower)
    m[:, :] = Diagonal(m)
    return
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

N = 512
SymmetricMatrix = rand(N, N)
SymmetricMatrix = SymmetricMatrix + SymmetricMatrix'
# SymmetricMatrix = RandTridiognalHermitian(N)
@time EigenDecomp = eigen(SymmetricMatrix)
@time RecursiveShiftedQR!(SymmetricMatrix)  # TridSymMatrix is destroyed because of mutability
display(sort(diag(SymmetricMatrix)))
display(EigenDecomp.values)
println("AbsoluteError: ", norm(EigenDecomp.values - sort(diag(SymmetricMatrix))))

# This is accurage, but the algorithm is not the impressive. :c