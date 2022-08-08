### 
### The inheritence graph is strictly a tree. Yeah. 

abstract type AbsSupSup1
end

abstract type AbsSupSup2
end

abstract type AbsSup <: Union{AbsSupSup1, AbsSupSup2}
end

mutable struct SupSup <: AbsSup
    function SupSup()
    return new() end
end

mutable struct SupSupU <: Union{AbsSupSup1, AbsSupSup2}
    function SupSupU()
    return new end
end

function Call_Me(this::AbsSup)
    typeof(this)|>display
return end

function Call_Me(this::SupSupU)
    typeof(this)|>display
return end

function Call_Me(this::AbsSupSup1)
    typeof(this)|>display
return end


