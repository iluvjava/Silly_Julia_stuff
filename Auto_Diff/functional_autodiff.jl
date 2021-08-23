### ============================================================================
### Differentialble Functions 
### ============================================================================

struct DiffUnivarFunc
    f::Function
    df::Function
    
    function DiffUnivarFunc(f, df)
        return new(f, df)
    end
end

function ∂(this::DiffUnivarFunc)
    return this.df
end

function Identity()
    return DiffUnivarFunc((x)-> x, (x) -> 1)
end

function Constant(c::Number)
    return DiffUnivarFunc((x)-> c, (x) -> 0)
end

function Base.:+(this::DiffUnivarFunc, other::DiffUnivarFunc)
    return DiffUnivarFunc(
        (x) -> this.f(x) + other.f(x),
        (x) -> this.df(x) + other.df(x)
    )
end

function Base.:-(this::DiffUnivarFunc)
    return DiffUnivarFunc(
        (x)-> - this.f(x), 
        (x)-> - other.df(x)
    )
end

function Base.:-(this::DiffUnivarFunc, that::DiffUnivarFunc)
    return DiffUnivarFunc(
        (x)-> - this.f(x) - that.f(x), 
        (x)-> - other.df(x) - that.df(x)
    )
end

function Base.:*(this::DiffUnivarFunc, that::DiffUnivarFunc)
    return DiffUnivarFunc(
        (x) -> this.f(x) + that.f(x),
        (x) -> begin x 
            this.f(x)*that.df(x) + 
            this.df(x)*that.f(x)
        end
    )
end

function Base.:*(this::Number, that::DiffUnivarFunc)
    return DiffUnivarFunc(
        (x) -> this*that.f(x),
        (x) -> this*that.df(x)
    )
end

function Base.:*(this::DiffUnivarFunc, that::Number)
    return that*this
end

function Base.:/(this::DiffUnivarFunc, that::DiffUnivarFunc)
    return DiffUnivarFunc(
        (x) -> this.f(x)/that.f(x), 
        (x) -> begin x
            -(this.f(x)'*that.df(x))/(that.f(x))^2 + 
            this.df(x)/that.f(x)
        end
    )
end

function Base.:∘(this::DiffUnivarFunc, that::DiffUnivarFunc)
    return DiffUnivarFunc(
        (x) -> (this.f∘that.f)(x), 
        (x) -> (this.df∘that.f)(x)*that.df(x)
    )
end

function Base.:^(this::DiffUnivarFunc, that::Number)
    return DiffUnivarFunc(
        (x) -> this.f(x)^that,
        (x) -> that*this.f(x)^(that - 1)*this.df(x)
    )
end

Sin = DiffUnivarFunc(
    sin, 
    cos
)

Cos = DiffUnivarFunc(
    cos, 
    (x) -> -sin(x)
)

Tan = DiffUnivarFunc(
    tan, 
    (x) -> sec(x)^2
)

Ln = DiffUnivarFunc(
    log, 
    (x) -> 1/x
)

Exp = DiffUnivarFunc(
    exp, 
    exp
)





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


