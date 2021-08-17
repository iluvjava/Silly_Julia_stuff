### **Intro**

We will be doing some examples Numerical Optimzation problems using Julia. 

**Packages**: 
* `JuMP`: [link](C:\Users\victo\source\repos\Silly_Julia_stuff\Mathematical Optimizations)
  * This one is more famous, and can interface with various different kind of solvers
  * This one has more tutorial examples. 
* `GLPK`: This is just a library in C and Julia Wraps around it. 
* `COSMO`: A pure julia optimizer, supports semi-definite programming, using the method of ADMM. 
* `JuliaOpt`: A collection of pure julia solver. Works with `convex.jl`, `JuMP.jl` interface. 
* `Convex.jl`: This is for convex problems interface, see [here](https://jump.dev/Convex.jl/stable/) for more. Usually people use stuff like: 
  * Second Order Cone representable objective 
    * SCS solver is the best
    * GLPK is for Linear programming, mixed integer programming. 
  * Linear objective
  * Advanced Features: 
    * Dual variables for each constraint. 
    * warmstart for minimization. 

---
### **A Simple Example**

Box projected L2 norm optimization. 

See [Simple_Example](./A%20Simple%20Example/Simple_Example.jl) for more. 

We use the 2-Norm matrix vector optimization problem using the `convex.jl`

The `convex.jl` interface is easier to manage and faster compare to the `JuMP.jl` interface. The later takes more time to construct a model. 

`JuMP.jl` is a lower level interface compare to the `Convex.jl` interface. For `JuMP.jl` we will need to write the loss function manually instead of using the builtin interface like `Convex.jl`

---
### **Lasso With L1 Norm**

Using the obe norm to do L1 regression analysis. 

> At the time of writing this, the `Convex.jl` interface doesn't support warmstarting for `SCS`, or the `COSMO` solver. For l1 lasso path, we write it in `JuMP.jl`, it's an lower 

Note: The loading of the `JuMP.jl` overhead might be outweighting the speed of the solvere. 


There are many other ways to solve the L1-Norm reglarization pronlem, But here we will be solving it using SCOP, we will treat it as a quadratic programming and then use ADMM based solver `COSMO.jl` .