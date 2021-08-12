### Use the lasso path for the crime train example. 

include("lasso.jl")

using DataFrames
using DelimitedFiles: readdlm

function DataSummary() 
    P, H = readdlm("./L1_Lasso_Module/crime-train.txt", '\t'; header=true)
    Df = DataFrame(P, reshape(H, length(H)))
    describe(Df)
end

println("")

DataSummary()




