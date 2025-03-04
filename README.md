double handed limit solver <br>
equation simplifier <br>
run the main.py file to start the application
```
>>> (1-cos(1-cos(x)))/(sin(x)^4)
((1+(-1*cos((1+(-1*cos(x))))))*(sin(x)^-4))
>>> limit 0
limit x->0 ((1+(-1*cos((1+(-1*cos(x))))))*(sin(x)^-4)) = (2^-3) = (2^-3)
>>> (sin(x)-tan(x))/x^3
((-1+cos(x))*(cos(x)^-1)*(x^-3)*sin(x))
>>> limit 0
limit x->0 ((-1+cos(x))*(cos(x)^-1)*(x^-3)*sin(x)) = (-2*(2^-1)*((2+(-1*(x^2)))^-1)) = (-2*(2^-2))
>>> (1-cos(100*x))/sin(x)^2
((1+(-1*cos((100*x))))*(sin(x)^-2))
>>> limit 0
limit x->0 ((1+(-1*cos((100*x))))*(sin(x)^-2)) = (10000*(2^-1)) = (10000*(2^-1))
>>> 
```
