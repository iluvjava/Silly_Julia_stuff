### ============================================================================
### Differential Number Pairs 
### This is more lightweight than function, but more copying if it's used too 
### intensly
### ============================================================================


struct DiffNum
    x::Number
    δx::Number

    function DiffNum(x, δ)
        return new(x, δ)
    end
    function DiffNum(x)
        return new(x, 1)
    end
end


function ∂(this::DiffNum)
    return this.δx
end


function Base.:+(this::DiffNum, that::DiffNum)
    return DiffNum(this.x + that.x, this.δx + this.δx)
end

function Base.:+(this::DiffNum, n::Number)
    return DiffNum(this.x + n, this.δx)
end

function Base.:+(n::Number, this::DiffNum)
    return DiffNum(this, n)
end


function Base.:*(this::DiffNum, that::DiffNum)
    return DiffNum(
            this.x*that.δx,
            this.δx*that.x + this.x*that.δx
        )
end

function Base.:*(this::DiffNum, n::Number)
    return DiffNum(
            n*this.x,
            n*this.δx
        )
end

function Base.:*(n::Number, this::DiffNum)
    return this*n
end


function Base.:/(this::DiffNum, that::DiffNum)
    return DiffNum(
        this.x/that.x, 
        -(this.x*that.δx)/that.x^2 + this.δx/that.x
    )
end

function Base.:/(this::DiffNum, n::Number)
    return DiffNum(
        this.x/n,
        this.δx/n
    )
end

function Base.:/(n::Number, this::DiffNum)
    return DiffNum(
        n/this.x,
        (-n/this.x^2)*this.δx
    )
end

function Base.:-(this::DiffNum, that::DiffNum)
    return DiffNum(
        this.x - that.x, 
        this.δx - that.δx
    )
end


function Base.:-(this::DiffNum)
    return DiffNum(
        -this.x,
        -that.δx
    )
end


function Base.:^(this::DiffNum, n::Number)
    return DiffNum(
        this.x^n, 
        n*this.x^(n - 1)
    )
end

function Base.:^(n::Number, this::DiffNum)
    return DiffNum(
        n^this.x, 
        n^this.x*log(n)
    )
end

function Base.sin(this::DiffNum)
    return DiffNum(sin(this.x), cos(this.x)*this.δx)
end

function Base.cos(this::DiffNum)
    return DiffNum(cos(this.x), -sin(this.x)*this.δx)
end

function Base.tan(this::DiffNum)
    return DiffNum(tan(this.x), this.δx*sec(this.x)^2)
end

function Base.exp(this::DiffNum)
    return DiffNum(exp(this.x), exp(this.x)*this.δx)
end

function Base.log(this::DiffNum)
    return DiffNum(log(this.x), this.δx/this.x)
end



### ----------------------------------------------------------------------------

### ============================================================================
### Differentiable multi-var func 
### Multiple inputs, one single output. 
### Multi-variable calculus Yeee ha 
### ============================================================================


struct DiffMultiNum
    x::Number
    ∇x::Vector{Number}
    
    function DiffMultiNum(x, ∇x)
        return new(x, ∇x)
    end

    function DiffMultiNum(x, dim::Int64)
        return new(x, ones(dim))
    end

end

