### **Intro**

We are interested in automatic differentiation in julia. 

> Composition of function value and their derivative should be enough. 

### **Functional Automatic Diff**

> Using lambda capture and higher order functions, we can make Automatic Differential Chain just by Lambda functions. 

* This is achievable by purly manipulating a pair of function and the derivative function. 
* And overriding numerical operators in the `Base` for our type: `DifferentiableFunction`. 