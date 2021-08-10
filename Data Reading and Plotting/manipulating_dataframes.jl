import DataFrames
import CSV

Df = CSV.read("programming_languages.csv", DataFrame)

display("The first row is ")
first(Df)

display("The last row is")
last(Df)

display("first 10 rows, of first columns, sliced out vector")
display(Df[1:10, "year"])

display("First 10 rows, of second column, slice out sub dataframe: ")
display(Df[1:10, ["language"]])

display("Slice out the second column as a vector")
display(Df[:, "language"])

display("another way of slicout out the vector")
display(Df.language)

display("Conditional row slice, language after the 1980: ")
display(Df[Df.year .> 1980, :])