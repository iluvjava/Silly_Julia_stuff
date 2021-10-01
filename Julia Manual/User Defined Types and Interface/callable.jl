mutable struct Callable
    field
    function Callable() 
        new(1)
    end
end

(::Callable)() = println("bruh, good call. ")
function (this::Callable)(arg::Any)
    println("bruh, good call with type $(typeof(arg)), this: $(this)")
end




