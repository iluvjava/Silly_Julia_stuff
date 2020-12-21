
println("mandebrot is running... ")
using Colors, Images, ImageView

function Manderbrot_Iterate(c0::Complex{Float64}; maxItr::Int64 = 1024, radius::Number = 2)
    zk = c0
    counter = 0
    while(abs(zk) < radius &&  counter < maxItr)
        zk = zk^2 + c0
        counter += 1
    end
    return counter # maxtr means in the set, else it's... not quiet in the set. 
end

"""
Given a region in the complex plane, discritized into a grids of a certain height, and width. 
We are interested in iterating all of them and make a coordinate iteration matrix 
out of it. 
"""
function Offset_Gridpoints(topLeft:: Tuple, bottomRight:: Tuple; 
                                width::Int64 = 1080, height::Int64 = 1920)
    
    dx = (bottomRight[1] - topLeft[1])/width
    dy = (topLeft[2] - bottomRight[2])/height
    Ygrid = topLeft[2] - dy/2: -dy: bottomRight[2] + dy/2
    Xgrid = topLeft[1] + dx/2: dx: bottomRight[1] - dx/2
    Xgrid = collect(Xgrid)
    Ygrid = collect(Ygrid)
    
    Ygrid = reshape(Ygrid, length(Ygrid), 1)
    Xgrid = reshape(Xgrid, 1, length(Xgrid))

    Ygrid = Ygrid .+ zeros(width, height)
    Xgrid = Xgrid .+ zeros(width, height)
    Plane = Xgrid + Ygrid.*im
    return Plane
end

ComplexPlane = Offset_Gridpoints((-2.5, 2), (1.5, -2); width=4000, height=4000)
ItrMatrix = Manderbrot_Iterate.(ComplexPlane; radius=4)
ItrMatrix = abs.(log.(ItrMatrix))
ItrMatrix = ItrMatrix./maximum(ItrMatrix)

img = Gray.(ItrMatrix)
save("Manderbrot.png", img)
println("mandebrot done.")

