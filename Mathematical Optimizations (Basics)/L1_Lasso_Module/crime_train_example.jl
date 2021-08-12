### Use the lasso path for the crime train example. 

include("lasso.jl")
include("utils.jl")

using DataFrames
using DelimitedFiles: readdlm
using Plots

function DataSummary() 
    P, H = readdlm("./L1_Lasso_Module/crime-train.txt", '\t'; header=true)
    Df = DataFrame(P, reshape(H, length(H)))
    return Df
end

PrintTitle("A Summarization of the Data. ")
TheData = DataSummary()
println(describe(TheData))
PrintTitle("We are gonna prepare the into a matrix and vector. ")
println("The number of predictors is: $(size(TheData, 1))")


function PrePareTheDataMatrix(data_frame::DataFrame) 
    """
        Min, max, median and mean are assumed to be important predictor
    """

end




