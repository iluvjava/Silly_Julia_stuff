module MyModule

    GlobalVar = 0

    function Func1()
        println("Function 1")
    end
    function Func2()
        println("Function 2")
        println("Global Variable is: $GlobalVar")
    end

    function ChangeModuleVarTo(n)
        global GlobalVar = n
    end

    export Func1, Func2, ChangeModuleVarTo
    export GlobalVar

    # The module will be added to main once it's evaluated. 
end
