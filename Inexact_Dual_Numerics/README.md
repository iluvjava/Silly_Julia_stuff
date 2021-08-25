### **Intro**

Trio Number computations can help us with inexact arithmetic, consider the following type of dual numbers: 


---
### **Error in Basic Arithematic**

To keep track of the error of a number, we would need both of these measures: 
* Relative Error 
* Absolute Error

> Multiplication will use the relative error and update the absolute error

> Addition and subtraction will use the absolute error and update the relative error. 

For any number, it has the following equivalent representation of errors involved: 

$$
a = 
\begin{cases}
    a \pm \delta 
    \\
    a(1 \pm \epsilon)
\end{cases} \text{where: } \epsilon = \frac{\delta}{|a|}
$$

read the code to see the update rule for of uncertainty for basic arithmetic. 


---
### **Orderness, Equality**

There are several way to do it: 

* Compare the $a$, original number and ignore the error bounds for the number. 
* Inpose strict comparison, give error if 2 numbers technically intersect with each other. 





