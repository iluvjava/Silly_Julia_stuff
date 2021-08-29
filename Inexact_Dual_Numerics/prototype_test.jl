using Test
include("prototype.jl")

function BasicSummation()

    for _ in 1:2^10
        inExactList = [InexactFloat(a*2.0^64) for a ∈ rand(2^4)]
        exactList = [Rational{BigInt}(a.x) for a ∈ inExactList]
        sum1 = sum(inExactList)
        sum2 = sum(exactList)
        # println("Inexact Number yeild: $(sum1)")
        # println("Exact Rational is: $(sum2)")
        # println("difference is: $(abs(sum1.x - float(sum2)))")
        # println("Is it in the error bound: $(float(sum2) ∈ sum1)")
        @assert float(sum2) ∈ sum1 "Ok one instance broke. "
    end
    return true
end

function BasicFiniteDiff()
    f(x) = 1 + 5x + 10x^2 + 10x^3 + 5x^4 + x^5

    Δx = eps(Float64)*2
    δf(x) = (f(x + Δx) - f(x))/Δx
    result1 = δf(InexactFloat(2.0))
    result2 = δf(Rational{BigInt}(2))
    @assert result2 ∈ result1 "Insnce Faile"* 
    "Number $(result2) is not in the error range of $(result1)"

    Δx = 1e-8
    δf(x) = (f(x + Δx) - f(x))/Δx
    result1 = δf(InexactFloat(2.0))
    result2 = δf(Rational{BigInt}(2))
    @assert result2 ∈ result1 "Insnce Faile"* 
    "Number $(result2) is not in the error range of $(result1)"

    return true

end

function BasicPolynomial()
    # Higher degree Polynomial
    p(x, deg=10) = sum([binomial(deg, KK)*x^KK for KK ∈ 0:deg]) # bad poly 
    q(x, deg=10) = (1 + x)^deg  # Good poly
    x = 1e-4
    @assert p(x) ∈ p(InexactFloat(x)) "This instance broke"
    @assert q(x) ∈ p(InexactFloat(x)) "This instance broke "
    return true
end

@testset "Basic" begin
    @test BasicSummation()
    @test BasicFiniteDiff()
    @test BasicPolynomial()
end


"""
    Does recursive summation asserts better numerical accuracy, according to
    our inexact number. 

"""
function TestRecursiveSummation()
    N = 2^22
    function RecursiveSum(
        arr::Vector, 
        ii::Union{Nothing, Int64}=nothing, 
        jj::Union{Nothing, Int64}=nothing
    )    
        if ii === nothing
            ii = 1
            jj = length(arr)
        end
        if jj - ii <= 1
            if ii == jj
                return arr[ii]
            end
            return arr[ii] + arr[jj]
        end

        m = (ii + jj)÷2
        leftSum = RecursiveSum(arr, ii, m)
        rightSum = RecursiveSum(arr, m + 1, jj)
        return leftSum + rightSum
    end
    toSum = [InexactFloat(II) for II in rand(N)]
    theSum1 = RecursiveSum(toSum)
    theSum2 = sum(toSum)
    theSum3 = sum([Rational{BigInt}(a.x) for a ∈ toSum])
    println("Using Recursive Sum: $(theSum1)")
    println("Using reduction Sum: $(theSum2)")
    println("Using the rational sum $(theSum3)")
    println("Abs and relative error for recursive and reduction sum are: ")
    for Sum ∈ [theSum1, theSum2]
        println("RelError: $((Sum.x - theSum3)/theSum3)")
        println("AbsError: $(Sum.x - theSum3)")
    end
    return true
end

@testset begin 

    @test TestRecursiveSummation()
end

 