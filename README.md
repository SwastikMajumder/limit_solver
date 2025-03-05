double handed limit solver <br>
equation simplifier <br>
run the main.py file to start the application
```
>>> 1+1
2
>>> sin(x)*sin(x)*sin(x)
(sin(x)^3)
>>> limit 0
limit x->0 (sin(x)^3) = 0
>>> (1-cos(x))/x^2
((1+(-1*cos(x)))*(x^-2))
>>> limit 0
limit x->0 ((1+(-1*cos(x)))*(x^-2)) = (2^-1)
>>> (sin(x)-tan(x))/x^3
((1+(-1*(cos(x)^-1)))*(x^-3)*sin(x))
>>> limit 0
limit x->0 ((1+(-1*(cos(x)^-1)))*(x^-3)*sin(x)) = (-3*(6^-1))
>>> sin(x)
sin(x)
>>> d/dx
cos(x)
>>> sin(sin(x)*cos(2*x))
sin((cos((2*x))*sin(x)))
>>> d/dx
(((-2*sin((2*x))*sin(x))+(cos((2*x))*cos(x)))*cos((cos((2*x))*sin(x))))
>>> (1-cos(100*x))/sin(x)^2
((1+(-1*cos((100*x))))*(sin(x)^-2))
>>> limit 0
limit x->0 ((1+(-1*cos((100*x))))*(sin(x)^-2)) = (10000*(2^-1))
>>> 
```
