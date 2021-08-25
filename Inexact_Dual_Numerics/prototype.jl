struct InexactFloat <: Real
    a::Float64 # Original number 
    δ::Float64 # Absolute error
    ϵ::Float64 # Relative Error

    function InexactFloat(a::Float64, error::Float64, mode::Symbol=:AbsoluteError)
        this = new()
        if mode == :AbsoluteError
            this.δ = abs(error)
            this.ϵ = abs(error)/abs(a)
            this.a = a
        elseif mode == :RelativeError
            this.ϵ = abs(error)
            this.δ = abs(error)*abs(a)
            this.a = a
        else
            error("Mode not recognized")
        end
        
        return this
    end
end

# Arithmetic operations for this given type ------------------------------------

function Base.:+(a::InexactFloat, b::InexactFloat)

end

function Base.:+(a::InexactFloat, b::Real)

end

function Base.:-(a::InexactFloat, b::InexactFloat)

end

function Base.:-(a::InexactFloat, b::Real)

end

function Base.:/(a::InexactFloat, b::InexactFloat) 

end

function Base.:/(a::InexactFloat, b::Real)

end

function Base.:/(a::Real, b::InexactFloat)

end

function Base.:^(a::InexactFloat, k::Real)

end


# Oderness, equality of the number ---------------------------------------------


# The special element in the field ---------------------------------------------

function Base.zero(a::InexactFloat)

end

function Base.show(io:IO, this::InexactFloat)
    
end