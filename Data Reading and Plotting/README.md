Reference Resources: [Here](https://github.com/JuliaAcademy/DataScience)

---
### **Intro**

In this environment, we will be learning about the basics of using Julia for data processing and visualization 

To read the data, we will need the following packages: 

* `DelimitedFiles.jl`: This is used for reading simple files with values that are separated by a character. [here](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/) for more. 
* `CSV.jl`: This is needed to read CSV files, and see [here](https://csv.juliadata.org/stable/) for more. 
* `XLSX.jl`: This is needed to parse excel files into data inside julia. See [here](https://felipenoris.github.io/XLSX.jl/stable/) for more. 

---
### **DataFrames**

The package that we are going to use is the [DataFrames.jl](https://dataframes.juliadata.org/stable/) package. 

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


---
### **Plotting**

You just need more practice to get use to it.

The **cheatsheet** is [here](./plotsjl-cheatsheet.pdf)

> However, do take notice that the plotting library for Julia has many different back ends. 

[Here](https://docs.juliaplots.org/latest/backends/) is a list of backend that is supported by `Plots` package in julia. 

And there is even the [Unicode Plot](https://github.com/Evizero/UnicodePlots.jl) which allows for plotting on the terminal. 
