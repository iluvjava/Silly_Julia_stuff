Reference Resources: [Here](https://github.com/JuliaAcademy/DataScience)

---
### **Intro**

In this environment, we will be learning about the basics of using Julia for data processing and visualization 


### **DataFrames**

reference [DataFrame Link](https://dataframes.juliadata.org/stable/)

Reading, writing the data: [here](./reading_and_writing.jl)

The manipulations of the dataframe is similar to Pandas in Python, but in Julia, we also have something more meta, that supports Linq like syntax for data frames manipulations, this is the: [DataFrameMeta](https://juliadata.github.io/DataFramesMeta.jl/stable/)

Basic Dataframe Manipulations: [here](./manipulating_dataframes.jl)


#### **DataFrames Manipulations: Split, Apply, Join**

This part is covered: [here](./data_manipulation_split_apply_combine.jl)

**Strategy: Gropuby and Combine** 

1. Specify the column you want to group that data frame by. 
   1. this will convert the dataframe to a Groupby Dataframe
   2. A groupy dataframe is a dataframe where one column is the same data for each dataframe. 
2. us the `combine` function to reduce on the gropuby data frame
   1. Each sub dataframe will be reduce to one row in the reduced dataframe
   2. We have the choice to do the following: 
        * Multiple columns to one: `mean`, `sum`, `average`
        * One columns to multiple: `minmax`, `min, max, std`
        * Multiple to Multiple: It's like a combinations of the 2 above. 

Use the do block it's the coolest and easiest to remember. 


