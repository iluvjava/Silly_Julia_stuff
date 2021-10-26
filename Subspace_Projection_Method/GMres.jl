include("iterative_hessenberg.jl")

# ==============================================================================
# Performs GMREs on Real Matrices
# ==============================================================================

mutable struct GMResReal
    ih::IterativeHessenberg
    R # Upper Tiangular Reduction form for the Hessenberg

    function GMResReal(A::Function, b, x0; max_k=typemax(Int64))
        this = new()
        this.ih = IterativeHessenberg(A, b; x0=x0)
        
        return this
    end

end

# ==============================================================================
# A class that performs plane rotations for the GRMes, for sake of speed, this
# won't actually use a matrix. 
# ==============================================================================

mutable struct CumulativePlaneRotation
    row_pairs::Vector
    initial_vector::Vector{Float64}

    function CumulativePlaneRotation()
        this = new()
        this.row_pairs = Vector{Vector{Reals}}()
        this.initial_vector = Vecotr{Reals}()
        return this
    end

end

"""
    Get the full rotation matrix for an instance of the 
    CumulativePlaneRotation 
"""
function GetFullRotationMatrix(this::CumulativePlaneRotation)
    
    
end

function GetFirstColumnOfRotationMatrix(this::CumulativePlaneRotation)
    n = length(this.row_pairs)
    Identity = Matrix{Float64}(I, n, n)
    # perform sussesive rotations on rows
    for IdxI ∈ 1: n - 1

    end
end

"""
    Load in sussesive 2 rows of the H matrix, and then store the 
    triangularization process for it. 
"""
function (this::CumulativePlaneRotation)(a::Real, b::Real)
    push!(this.row_pairs, Vector{reals}())
    append!(this.row_pairs[end], [a, b])
    return nothing
end



