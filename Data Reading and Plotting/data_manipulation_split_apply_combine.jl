using DataFrames, CSV, Statistics

Block1 = begin
    Iris = CSV.read("iris.csv", DataFrame)
    GroupBySpecies = groupby(Iris, :Species)
end

# Reduction on one attribute to one attribute
Block2 = begin 
                           # data frame, Keyitem 
    CombineToMean = combine(GroupBySpecies, :PetalLength => mean)
end

# Reduction on multiple attributes to multiple atttribute
Block3 = begin
    CombineMultiToMulti = combine(GroupBySpecies, 
        [:PetalLength, :SepalLength] => ((p, s) -> (MeanPetallength=mean(p)/mean(s), SepalLengthSum=sum(p))) => AsTable)
        # [:col1, :col2] => ((col1, col2) -> (func(col1, col2), func(col1, col2)) => As table
        # It's a nested keyitem value pair, this is something programmer should remember. 
end

# Reduction on multiple attributes to one attributes 
# x: each sub dataframe 
Block3 = begin 
    CombineMultiToOne = combine(x -> std(x.PetalLength) / std(x.SepalLength), GroupBySpecies)
    # x::DataFrame -> ::ColumnType, Gropubed by data
end


# Reduction transforming one attribute to multiple
Block4 = begin
    combine(GroupBySpecies,:PetalLength => (x -> [extrema(x)]) => [:min, :max])
           #GropubyTable, :col1 => (col1 -> [Tuple]) ==> [func1, func2]
end 

# Reduction on gropuby one attribute by one attribute
Block5 = begin 
    combine(GroupBySpecies, nrow, :PetalLength => mean => :mean)
end


# Finally let's check out a simpler syntax for reduction involving multiple 
# columns on the gropuby function

Block6 = begin 
    combine(GroupBySpecies) do df # alis for the Sub-Dataframe
        (m = mean(df.PetalLength), s² = var(df.PetalLength))
    end
end