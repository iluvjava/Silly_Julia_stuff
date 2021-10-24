include("iterative_hessenberg.jl")

# ==============================================================================
# Performs GMREs on Real Matrices
# ==============================================================================

mutable struct GMResReal
    ih::IterativeHessenberg
    R # Upper Tiangular Reduction form for the Hessenberg



end

# ==============================================================================
# A class that performs plane rotations for the GRMes, for sake of speed, this
# won't actually use a matrix. 
# ==============================================================================

mutable struct CumulativePlaneRotation
    row_pairs::Vector

    function CumulativePlaneRotation()
        this = new()
        this.row_pairs = Vector{Vector{Reals}}()
        return this
    end

end

"""
    Get the full rotation matrix for an instance of the 
    CumulativePlaneRotation 
"""
function GetFullRotationMatrix(this::CumulativePlaneRotation)
    
    
end



