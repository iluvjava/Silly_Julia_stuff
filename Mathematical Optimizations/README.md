### **Intro**

We will be doing some examples Numerical Optimzation problems using Julia. 

**Packages**: 
* `JuMP`: [link](C:\Users\victo\source\repos\Silly_Julia_stuff\Mathematical Optimizations)
  * This one is more famous, and can interface with various different kind of solvers
  * This one has more tutorial examples. 
* `GLPK`: This is just a library in C and Julia Wraps around it. 
* `Convex.jl`: This is for convex problems, see [here](https://jump.dev/Convex.jl/stable/) for more. Usually people use stuff like: 
  * Second Order Cone representable objective 
    * SCS solver is the best
    * GLPK is for Linear programming, mixed integer programming. 
  * Linear objective
  * Advanced Features: 
    * Dual variables for each constraint. 
    * warmstart for minimization. 

---
### **A Simple Example**


