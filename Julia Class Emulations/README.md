### **Intro**

In this file, we are going to make some java-like classes in Julia. 

This is possible via lamabda closure in Julia. 


---
### **Explanation** 

* Treat of function is object like in julia, it can be defined inside a function. 

*The inner function capture the variabless in the outter scope. 

* Using this, we can make function fields for a given type, and then assign the function in the constructor, which is possible via `new()`, a function that returns the current instance of the object. 

* And we will need to explicity return the new instance of the object in the constructor, or else we lose the reference to the object. 


Pretty cool huh? 
