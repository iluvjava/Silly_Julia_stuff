include("core_optimized.jl")

using BenchmarkTools

function TestReconstructionError()
    A = rand(5, 5) .- 0.5
    A += rand(5, 5)*im
    H = similar(A, ComplexF64)
    Q = similar(A, ComplexF64)
    ArnoldiIterate!(A, Q, H)
    println("Q Should be Orthogonal")
    display(Q*Q')
    println("QHQ' should be similar to A")
    println("This is Q*H*Q'")
    display(Q*H*Q')
    println("This is A")
    display(A)
    println("Reconstruction Error: ")
    Error = sum(abs.(Q*H*Q' - A))
    println(Error)
    @assert Error <= 1e-6 "error should be small"
    print("Test Ok")
end

function SpeedTest()
    N = 300
    rep = 10
    A = rand(N, N) .- 0.5
    A += (rand(N, N) .- 0.5).*im
    H = similar(A, ComplexF64)
    Q = similar(A, ComplexF64)  
    timetook = @elapsed for _ in 1:rep   
        Q, H = ArnoldiIterate!(A, Q, H)
    end
    println(string("Total time took to decompose ",
                    rep, 
                    ", " , N, " by ", N, " matrices",
                    " is: ", timetook))
    println("This is for the optimized version of the code. ")
end


TestReconstructionError()
SpeedTest()