include("mytypes.jl")
function test()
    TestType = Int128
    tree = BinaryTree{TestType}()

    @time for II in rand(1:typemax(TestType), 2^20)
        push!(tree, convert(TestType, II))
    end

    PreviousItem = -Inf
    for Item in tree.root
        if Item < PreviousItem
            println(string(item, " < ", PreviousItem))
            throw(AssertionError("Something is wrong. "))
        end
        PreviousItem = Item
    end
    println("A simple Test past. ")
end

test() 