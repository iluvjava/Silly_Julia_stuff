# This file is about basic data manipulations

using BenchmarkTools
using DataFrames      # for construction of dataframe. 
using DelimitedFiles  # For delimited files
using CSV             # for CSV
using XLSX            # for Excel Spread sheet
using Downloads       # for fetch from URL 

# Download into pwd and return a string of the file name
P = Downloads.download("https://raw.githubusercontent.com/nassarhuda"*
                        "/easy_data/master/programming_languages.csv", 
                       "programming_languages.csv")


# ------------------------------------------------------------------------------
# read a file using the DelimitedFiles Library 
# Header and body is read separately, with delimiter set to ",". 
P, H = readdlm("programming_languages.csv", ','; header=true);

# file name, data, delimiter, and options
writedlm("programming_languages.txt", P, "-")


# Read if we already know that this is a CSV file: 
# file name, datatype
C = CSV.read("programming_languages.csv", DataFrame);
println("DataFrame: ")
display(C)
println("Header Names")
println(names(C))
println("Describe()")
display(describe(C))

# Storing to a CSV file, with a matrix P 
# file name, an instance of dataframe. 
# Use ?CSV.write on julia REPL for more information. 
CSV.write("programminglange_CSV.csv", DataFrame(P, :auto))