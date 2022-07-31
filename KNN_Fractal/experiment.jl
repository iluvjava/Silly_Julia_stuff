struct CantorDust
    level::Int
    cut_out::Number

    function CantorDust(level=4, cut_out=1/3)
        if cut_out > 1 || cut_out < 0
            error("Error, invalid cut_out parameter, it should be between 0 and 1. ") 
        end
        this = new(level, cut_out) 
    return this end
end

function Iterate(this::CantorDust, p::Number)
    if p < 0 || p > 1
        return false
    end
    l = this.level
    c = this.cut_out
    while l > 0
        if p > c && p < 1 - c
            return this.level - l
        end
        if p <= c
            p *= 1/c
        else
            p = 1 + (p - 1)*(1/c)
        end
        l -= 1
    end
return this.level end


function (this::CantorDust)(a::Number, b::Number)
return min(Iterate(this, a), Iterate(this, b)) end



### Experiments
using Plots, Images, ImageShow

N = 3^8
l = 30
C = CantorDust(l)
step = 1/N
start = 1/(2*N)
final = 1 - step
grid = [C(x, y) for x in start:step:final, y in start:step:final] .|> Float64
grid = grid/maximum(grid) .- minimum(grid)
grid = (grid .+ 1) .|> log2
# "dust count: $(sum(grid)/N^2)"|>println
save("interesting_dust.png", grid |> colorview(Gray))

"Done" |> println

