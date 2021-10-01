mutable struct Callable


end

(::Callable)() = println("bruh, good call. ")
function (::Callable)(arg::Any)
    println("bruh, good call with type $(typeof(arg))")
end

function (::Callable)func(arg::Any)
    println("bruh, good call with type $(typeof(arg))")
end

