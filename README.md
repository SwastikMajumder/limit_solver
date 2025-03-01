double handed limit solver <br>
equation simplifier <br>
run the main.py file to start the application
```
>>> x+x+x*x
((2*x)+(x^2))
>>> (sin(x)*sin(x))/sin(x)+1+1+1+1+1*x*x*x*x+(sin(x)*sin(x))/sin(x)+1+1+1+1+1*x*x*x*x+(sin(x)*sin(x))/sin(x)+1+1+1+1+1*x*x*x*x
(12+(3*(x^4))+(3*sin(x)))
>>> sin(x*y)/sin(x*z)
((sin((x*z))^-1)*sin((x*y)))
>>> limit 0
limit x->0 ((sin((x*z))^-1)*sin((x*y))) = ((z^-1)*y)
>>> x/sin(10*x)
((sin((10*x))^-1)*x)
>>> limit 0
limit x->0 ((sin((10*x))^-1)*x) = (10^-1)
>>> x/cos(x)
((cos(x)^-1)*x)
>>> limit 0
limit x->0 ((cos(x)^-1)*x) = 0
>>>
```
