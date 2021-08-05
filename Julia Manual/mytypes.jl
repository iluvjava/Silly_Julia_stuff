### ----------------------------------------------------------------------------
###  This is a basic Binary Search Tree Class with some methods
### ----------------------------------------------------------------------------

using BenchmarkTools
using InteractiveUtils


writeLine(x...)=(println ∘ string)(x...)


# An abstract type, that is like, purely a type and no other usage, you can't
# instantiante it. 
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
        new{T}(nothing, nothing, nothing, item)
    end
end


# ----------------------- Binary Search tree -----------------------------------

mutable struct BinaryTree{T} <: SortedSet{T}
    # It has a root to it
    root::Union{BinaryNode{T}, Nothing}
    count::UInt32
    
    function BinaryTree{T}() where T
        new{T}(nothing, 0)
    end
end


function Base.push!(this::BinaryTree{T}, item::T)::Nothing where {T}
    Added::Bool = false
        function insert!(
                this::Union{BinaryNode{T}, Nothing},
                parent::Union{BinaryNode{T}, Nothing},
                item::T
            )::BinaryNode{T} where T
            
            if this === nothing
                NewNode = BinaryNode{T}(item)
                NewNode.parent = parent
                Added = true # modifies the outter scope 
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
# States: 
#   Current Node
#   Integer of {1, 2, 3}
#       0. printing left 
#       1. printing middle
#       2. printing right
#       3. ignore the node

function Base.iterate(this::BinaryNode)
    return iterate(this, (this, [0x0]))
end

function Base.iterate(this::BinaryNode, states_chain::Tuple)
    # Assign types to speed it up. 
    n = states_chain[1]::BinaryNode
    v = states_chain[2]::Vector{UInt8}
    while length(v) != 0
        # root node
        if length(v) == 0
            return (this.data, states_chain)
        end
        s = v[end]
        if s == 0
            v[end] = 0x1
            if n.left === nothing
                continue
            else
                append!(v, 0)
                n = n.left
            end
        elseif s == 1
            v[end] = 0x2
            return (n.data, (n, v))
        elseif s == 2
            v[end] = 0x3
            if n.right === nothing 
                continue
            else
                append!(v, 0)
                n = n.right
            end
        else
            pop!(v)
            n = n.parent
        end
    end
    return nothing
end

