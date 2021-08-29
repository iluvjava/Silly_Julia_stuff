function QuantizeError(a::Float64)
    return 2.0^(exponent(a) - 54)
end