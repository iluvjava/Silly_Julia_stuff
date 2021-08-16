# We can use the julia language to reproduce OOP-like
# bahaviors using the lambda capturing functionality. 
# 
# Ok, let's us this example to make a point here. 


mutable struct MyCartesianPoint
    x::Float64
    y::Float64

    GetX::Function
    GetY::Function
    SetX::Function
    SetY::Function 
    function MyCartesianPoint(x::T, y::S) where {T<:Real, S<:Real}
        this = new()  # The instance
        this.x = x
        this.y = y

        SetX(x::T) = begin
            this.x = x
            return nothing
        end
        SetY(y::S) = begin
            this.y = y
            return nothing
        end 
        this.SetX = SetX
        this.SetY = SetY
        this.GetX = () -> this.x
        this.GetY = () -> this.y
        
        return this  # You actually need to do this explcitly
    end

end

# And, a way of binding that closure with static runctions are also a pretty good approach. 




