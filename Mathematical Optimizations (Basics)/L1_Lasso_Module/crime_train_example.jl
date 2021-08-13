### Use the lasso path for the crime train example. 

include("lasso.jl")
include("utils.jl")

using DataFrames
using DelimitedFiles: readdlm
using Plots
import UnicodePlots as Uplot


function DataSummary() 
    P, H = readdlm("./L1_Lasso_Module/crime-train.txt", '\t'; header=true)
    Df = DataFrame(P, reshape(H, length(H)))
    return Df
end


PrintTitle("A Summarization of the Data: The Predictors")
TheData= DataSummary()
println(describe(TheData))
println("The size of the dataframe is: $(size(TheData))")
println("The predictant is:\"$(names(TheData)[1])\"")
PrintTitle("We are gonna prepare the into a matrix and vector. ")

function PrepareTheDataMatrix(data_frame::DataFrame) 
    """
        Min, max, median and mean are assumed to be important predictor
    """
    FeatureMatrix = data_frame[!, 2:end]
    FeatureMatrix = Matrix(data_frame)
    Predictants = data_frame[!, 1]
    return FeatureMatrix, Predictants
end

A, b = PrepareTheDataMatrix(TheData)

function AnalyaisWithLasso(A::Matrix, b::Matrix )

