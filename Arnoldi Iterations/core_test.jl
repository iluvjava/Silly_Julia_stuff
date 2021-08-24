include("core.jl")

using BenchmarkTools

function SimpleTest()
    A = rand(5, 5) .- 0.5
    A += rand(5, 5)*im
    Q, H = ArnoldiIterate(A)
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
    @time timetook = @elapsed for _ in 1:rep   
        Q, H = ArnoldiIterate(A)
    end
    println(string("Total time took to decompose ",
                    rep, 
                    ", " , N, " by ", N, " matrices",
                    " is: ", timetook))
    println("This is for the un-optimized version of the code. ")
end

SimpleTest()
SpeedTest()

