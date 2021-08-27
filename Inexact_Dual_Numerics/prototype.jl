struct InexactFloat <: Real
    x::Float64 # Original number 
    δ::Float64 # Absolute error
    ϵ::Float64 # Relative Error

    function InexactFloat(a::Float64)
        ϵ = 2^(-54)
        δ = abs(a)*ϵ
        x = a
        return new(x, δ, ϵ)
    end
    
    function InexactFloat(a::Float64, ϵ::Float64)
        x = abs(a)
        ϵ = max(abs(ϵ), 2^(-54))
        δ = ϵ*abs(a)
        return new(x, δ, ϵ)
    end
end


# Arithmetic operations for this given type ------------------------------------

function Base.:+(a::InexactFloat, b::InexactFloat)
    if a.x < b.x
        return b + a  # WLOG: a > b
    end
    c = a.x + b.x
    δ = a.δ + b.δ + (b.x - (c - a.x)) # compensation, kahan summation
    ϵ = δ/abs(c)
    return InexactFloat(c, ϵ)
end
function Base.:+(a::InexactFloat, b::Real)
    c = a.x + b
    δ = a.δ
    ϵ = δ/abs(c)
    return InexactFloat(c, ϵ)
end
function Base.:+(b::Real, a::InexactFloat)
    c = a.x + b
    δ = a.δ
    ϵ = δ/abs(c)
    return InexactFloat(c, ϵ)
end

function Base.:-(a::InexactFloat, b::InexactFloat)
    c = a.x - b.x
    δ = a.δ + b.δ
    ϵ = δ/abs(c)
    return InexactFloat(c, ϵ)
end
function Base.:-(a::InexactFloat, b::Real)
    c = a.x - b
    δ = a.δ
    ϵ = δ/abs(c)
    return InexactFloat(c, ϵ)
end
function Base.:-(b::Real, a::InexactFloat)
    c = a.x - b
    δ = a.δ
    ϵ = δ/abs(c)
    return InexactFloat(c, ϵ)
end
function Base.:-(a::InexactFloat)
    c = -a.x 
    ϵ = a.ϵ
    return InexactFloat(c, ϵ)
end

function Base.:*(a::InexactFloat, b::InexactFloat)
    c = a.x*b.x
    ϵ = a.ϵ + b.ϵ + a.ϵ*b.ϵ
    return InexactFloat(c, ϵ)
end
function Base.:*(a::InexactFloat, b::Real)
    c = a.x*b
    ϵ = a.ϵ
    return InexactFloat(c, ϵ)
end
function Base.:*(b::Real, a::InexactFloat)
    c = a.x*b
    ϵ = a.ϵ
    return InexactFloat(c, ϵ)
end

function Base.:/(a::InexactFloat, b::InexactFloat) 
    c = a.x/b.x
    if b.ϵ >= 1
        error("Inexact float is probably dividing by zero. ")
    end
    ϵ = abs((a.ϵ + b.ϵ)/(1 - b.ϵ))  # Catastropic cancellation
    return InexactFloat(c, ϵ)
end
function Base.:/(a::InexactFloat, b::Real)
    c = a.x/b
    ϵ = a.ϵ
    return InexactFloat(c, ϵ)
end
function Base.:/(a::Real, b::InexactFloat)
    c = a.x + b
    ϵ = b.ϵ/(1 + b.ϵ)
    return InexactFloat(c, ϵ)
end

function Base.:^(a::InexactFloat, k::Integer)
    # Big error no worry about cancellation. 
    if a.ϵ <= 2^(-8)
        return InexactFloat(a.x^k, (1 + a.ϵ)^k - 1)
    end
    # small error we need to be smarter. 
    ϵ = exp(k*a.ϵ) - 1
    c = a.x^k
    return InexactFloat(c, ϵ)
end

# Oderness, equality of the number ---------------------------------------------


# Comon Functions on this type of numbers --------------------------------------


# The special element in the field ---------------------------------------------

function Base.zero(a::InexactFloat)
    return InexactFloat(0, 0)
end

function Base.show(io::IO, this::InexactFloat)
    show(io, "$(this.x) ± $(this.δ)")
end