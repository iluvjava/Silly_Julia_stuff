### https://discourse.julialang.org/t/redefine-struct-when-working-with-repl/25942

module M
    mutable struct MyStruct
        attribute1::Int
        attribute2::Int
        # attribute3::Int # Change attribute here and loading it won't end the REPL Sessions. 
        function MyStruct() 
        return new() end
    end
    export MyStruct
end



