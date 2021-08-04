using BenchmarkTools
using InteractiveUtils


writeLine(x...)=(println ∘ string)(x...)


# An abstract type, that is like, purely a type and no other usage, you can't instantiante it. 
abstract type SortedSet{T}

end

# ---------------------------- binary Node -------------------------------------

mutable struct BinaryNode{T}

    # for iterator 
    parent::Union{BinaryNode{T}, Nothing}

    # Nullable type
    left::Union{BinaryNode{T}, Nothing}
    right::Union{BinaryNode{T}, Nothing}
    data:: T # non nullable type

    # This disable the default constructor. 
    function BinaryNode{T}(item) where T 
        new(nothing, nothing, nothing, item)
    end
end


# ----------------------- Binary Search tree -----------------------------------
mutable struct BinaryTree{T} <: SortedSet{T}
    # It has a root to it
    root::Union{BinaryNode{T}, Nothing}
    count::UInt32
    
    function BinaryTree{T}() where T
        new(nothing, 0)
    end
end


function Base.push!(this::BinaryTree{T}, item::T)::Nothing where T
    Added::Bool = false
        function insert!(
                this::Union{BinaryNode{T}, Nothing},
                parent::Union{BinaryNode{T}, Nothing},
                item::T
            )::BinaryNode{T} where T
            
            if this === nothing
                NewNode = BinaryNode{T}(item)
                NewNode.parent = parent
                Added = true     
                return NewNode
            end

            if this.data === item  
                # don't put this in. 
                return this
            end

            if item <= this.data
                this.left = insert!(this.left, this, item)
            else
                this.right = insert!(this.right, this, item)
            end

            return this
        end
    
    
    this.root = insert!(this.root, nothing ,item)      
    if Added 
        this.count += 1
    end
    return nothing # has to be here or it will return to return this.root
end

function Base.length(this::BinaryTree)::UInt32
    return this.count
end

# ------------------------- Iterate --------------------------------------------
# The traversal will be in order 
function Base.iterate(this::BinaryTree{T}) where T
    if this.root === nothing
        return nothing
    end
    # start with the smallest 
    while !(this.left === nothing)
        this = this.left
    end
    return (this.data, this)
end

function Base.iterate(A::BinaryTree{T}, state::BinaryNode{T}) where T
    


end

tree = BinaryTree{Int64}()
push!(tree, 3)
push!(tree, 4)
push!(tree, 2)
println(tree.root.data)
println(tree.root.left.data)
println(tree.root.right.data)
println(length(tree))

