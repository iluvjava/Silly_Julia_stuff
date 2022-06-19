struct CantorDust
    level::Int

    function CantorDust(level=4)
        this = new(level) 
    return this end

    
end

function InCantor(this::CantorDust, p::Number)
    if p < 0 || p > 1
        return false
    end

    l = this.level
    while l >= 0
        if p > 1/3 && p < 2/3
            return false
        end
        if p <= 1/3
            p *= 3
        else
            p = 1 + (p - 1)*3
        end
        l -= 1
    end
return true end


function (this::CantorDust)(a::Number, b::Number)
return InCantor(this, a)&&InCantor(this, b) end



### Experiments
using Plots, Images, ImageShow

N = 4096

for l in 2:10
    C = CantorDust(l)
    grid = [C(x, y) for x in LinRange(0, 1, N), y in LinRange(0, 1, N)].|>Float64
    "dust count: $(sum(grid)/N^2)"|>println
    save("dust_$(l).png", colorview(Gray, grid))
end


