### The absolute error is quantized. 

struct InexactFloat <: AbstractFloat
    x::Float64    # Original number 
    δ::Float64    # Absolute error
    ϵ::Float64    # Relative Error

    function InexactFloat()
        return InexactFloat(0.0)
    end
    
    function InexactFloat(a::Float64)
        if a == 0.0
            δ = 0
        else
            δ = 2.0^(exponent(a) + exponent(eps(Float64)) - 1)
        end
        ϵ = δ/abs(a)
        x = a
        return new(x, δ, ϵ)
    end
    
    function InexactFloat(a::Float64, ϵ::Float64)
        x = a
        ϵ = max(abs(ϵ), eps(Float64)/2)
        δ = abs(a)*ϵ
        return new(x, δ, ϵ)
    end
end

# Arithmetic operations for this given type ------------------------------------

function Base.:+(a::InexactFloat, b::InexactFloat)
    if a.x < b.x
        return b + a  # WLOG: a > b
    end
    c = a.x + b.x
    # Loss of significance cancellation.
    δ = a.δ + b.δ + 2.0^(exponent(a.x) + exponent(eps(Float64)) - 1)
    ϵ = δ/abs(c)
    return InexactFloat(c, ϵ)
end


function Base.:+(a::InexactFloat, b::Real)
    c = a.x + b
    return InexactFloat(c, a.ϵ)
end


function Base.:+(b::Real, a::InexactFloat)
    c = a.x + b
    return InexactFloat(c, a.ϵ)
end


function Base.:-(a::InexactFloat, b::InexactFloat)
    c = a.x - b.x
    δ = a.δ + b.δ + 2.0^(exponent(a.x) + exponent(eps(Float64)) - 1)
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
    ϵ = a.ϵ + 2.0^(exponent(eps(Float64)) - 1) + 
        a.ϵ*2.0^(exponent(eps(Float64)) - 1)
    return InexactFloat(c, ϵ)
end


function Base.:*(b::Real, a::InexactFloat)
    c = a.x*b
    ϵ = a.ϵ + 2.0^(exponent(eps(Float64)) - 1) + 
    a.ϵ*2.0^(exponent(eps(Float64)) - 1)
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

# TODO: Integer Division
# TODO: Power is also Inexact 


# Oderness, equality of the number ---------------------------------------------
# This will be needed to run the linear algebra libraries and other common
# Numerical Algorithms. 


# Comon Functions on this type of numbers --------------------------------------




# The special functons ---------------------------------------------------------

function Base.zero(::Type{InexactFloat})
    return InexactFloat(0.0, 0.0)
end

function Base.show(io::IO, this::InexactFloat)
    show(io, "$(this.x) ± $(this.δ)")
end

function Base.:∈(a::Real, b::InexactFloat)
    return a <= b.x + b.δ && a >= b.x - b.δ
end

function Base.:∈(a::InexactFloat, b::InexactFloat)
    return a.x + a.δ <= b.x + b.δ && a.x - a.δ >= b.x - b.δ
end

function Base.:≈(a::Real, b::InexactFloat)
    return a ∈ b
end

function Base.convert(::Type{InexactFloat}, a::Float64)
    return InexactFloat(a)
end
function Base.convert(::Type{InexactFloat}, a::InexactFloat)
    return a
end
function Base.convert(::Type{InexactFloat}, a::Integer)
    return InexactFloat(a + 0.0)
end

function Base.convert(::Type{Float64}, a::InexactFloat)
    return a.x
end
function Base.convert(::Type{Integer}, a::Integer)
    return convert(Integer, a.x)
end

