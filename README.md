symbolic mathematics solver <br>
no libraries used for building this algorithm (except for using lark parser for parsing math equation into expression trees) <br>
run the main.py file to start the application
```
>>> (x-1)/(x+3)
((-1+x)*((3+x)^-1))
>>> d/dx
((-1*(-1+x)*((3+x)^-2))+((3+x)^-1))
>>> fraction
((12+(4*x))*((3+x)^-3))
>>> factor
(4*((3+x)^-2))
>>> sin(x^2)*x
(sin((x^2))*x)
>>> Sdx
((-2^-1)*cos((x^2)))
>>> d/dx
(sin((x^2))*x)
>>> sin(cos(x^6+1))
sin(cos((1+(x^6))))
>>> d/dx
(-6*cos(cos((1+(x^6))))*(x^5)*sin((1+(x^6))))
>>> (sin(x)-tan(x))/x^3
(((-1*(cos(x)^-1)*sin(x))+sin(x))*(x^-3))
>>> limit 0
limit x->0 (((-1*(cos(x)^-1)*sin(x))+sin(x))*(x^-3)) = (-2^-1)
>>> sin(2*x)/x
((x^-1)*sin((2*x)))
>>> limit 0
limit x->0 ((x^-1)*sin((2*x))) = 2
>>> 1/tan(x)-tan(x)-2/tan(2*x)
((-1*(cos(x)^-1)*sin(x))+(-2*((cos((2*x))^-1)^-1)*(sin((2*x))^-1))+(((cos(x)^-1)^-1)*(sin(x)^-1)))
>>> trig
((-1*((-1*(sin(x)^2))+(cos(x)^2))*(cos(x)^-1)*(sin(x)^-1))+(-1*(cos(x)^-1)*sin(x))+(cos(x)*(sin(x)^-1)))
>>> fraction
0
>>> (1-cos(2*x))/sin(2*x)-tan(x)
((-1*(cos(x)^-1)*sin(x))+((1+(-1*cos((2*x))))*(sin((2*x))^-1)))
>>> trig
((-1*(cos(x)^-1)*sin(x))+((2^-1)*(1+(-1*((-1*(sin(x)^2))+(cos(x)^2))))*(cos(x)^-1)*(sin(x)^-1)))
>>> fraction
((2^-1)*(cos(x)+(-1*(cos(x)+(-1*cos(x)*(sin(x)^2))))+(-1*cos(x)*(sin(x)^2)))*((cos(x)^2)^-1)*(sin(x)^-1))
>>> expand
0
>>> (1+cos(2*x))/sin(2*x)-1/tan(x)
((-1*((cos(x)^-1)^-1)*(sin(x)^-1))+((1+cos((2*x)))*(sin((2*x))^-1)))
>>> trig
((-1*cos(x)*(sin(x)^-1))+((2^-1)*(1+(-1*(sin(x)^2))+(cos(x)^2))*(cos(x)^-1)*(sin(x)^-1)))
>>> fraction
((2^-1)*((-1*(1+(-1*(sin(x)^2)))*sin(x))+(-1*(sin(x)^3))+sin(x))*(cos(x)^-1)*((sin(x)^2)^-1))
>>> fraction
((2^-1)*((-1*(1+(-1*(sin(x)^2)))*sin(x))+(-1*(sin(x)^3))+sin(x))*(cos(x)^-1)*(sin(x)^-2))
>>> expand
0
>>> sin(2*x)/(1+cos(2*x))-tan(x)
((-1*(cos(x)^-1)*sin(x))+(((1+cos((2*x)))^-1)*sin((2*x))))
>>> trig
((-1*(cos(x)^-1)*sin(x))+(2*cos(x)*((1+(-1*(sin(x)^2))+(cos(x)^2))^-1)*sin(x)))
>>> expand
((-1*(cos(x)^-1)*sin(x))+(2*cos(x)*((1+(-1*(sin(x)^2))+(cos(x)^2))^-1)*sin(x)))
>>> fraction
(((-1*sin(x))+((1+(-1*(sin(x)^2)))*sin(x))+(sin(x)^3))*((2+(-2*(sin(x)^2)))^-1)*(cos(x)^-1))
>>> expand
0
>>> (x+5)*(2*x+6)
((5+x)*(6+(2*x)))
>>> expand
(30+(16*x)+(2*(x^2)))
>>> factor
(2*(3+x)*(5+x))
>>> tan(3*x)*tan(5*x)*tan(8*x)
((cos((3*x))^-1)*(cos((5*x))^-1)*(cos((8*x))^-1)*sin((3*x))*sin((5*x))*sin((8*x)))
>>> Sdx
(((-8^-1)*log(abs(cos((8*x)))))+((3^-1)*log(abs(cos((3*x)))))+((5^-1)*log(abs(cos((5*x))))))
>>> tan(x*6+7)
((cos((7+(6*x)))^-1)*sin((7+(6*x))))
>>> Sdx
((-6^-1)*log(abs(cos((7+(6*x))))))
>>> sin(x)/sin(x-z)
((sin(((-1*z)+x))^-1)*sin(x))
>>> Sdx
((((-1*z)+x)*cos(z))+(log(abs(sin(((-1*z)+x))))*sin(z)))
>>> expand
((-1*cos(z)*z)+(cos(z)*x)+(log(abs(sin(((-1*z)+x))))*sin(z)))
>>> sin(x+x+x+x+x*x+x^2+x+x+x+x)*sin(2+x)*sin(1+x+1)
((sin((2+x))^2)*sin(((2*(x^2))+(8*x))))
>>> d/dx
((2*cos((2+x))*sin((2+x))*sin(((2*(x^2))+(8*x))))+((8+(4*x))*cos(((2*(x^2))+(8*x)))*(sin((2+x))^2)))
>>> 1
1
>>> 1+4+5+6
16
>>> d/dx
0
>>> 1/sin(x)^2+1/cos(x)^2
((cos(x)^-2)+(sin(x)^-2))
>>> trig
((cos(x)^-2)+(sin(x)^-2))
>>> fraction
(((1+(-1*(sin(x)^2)))^-1)*((sin(x)^2)^-1))
>>> fraction
(((1+(-1*(sin(x)^2)))^-1)*(sin(x)^-2))
>>> expand
(((1+(-1*(sin(x)^2)))^-1)*(sin(x)^-2))
>>> (x+1)*(x-1)
((-1+x)*(1+x))
>>> expand
(-1+(x^2))
>>> factor
((-1+x)*(1+x))
>>> (1/sin(x)-1/tan(x))^2-(1-cos(x))/(1+cos(x))
((-1*(1+(-1*cos(x)))*((1+cos(x))^-1))+(((-1*((cos(x)^-1)^-1)*(sin(x)^-1))+(sin(x)^-1))^2))
>>> fraction
((-1+cos(x)+(((-1*(1+(-1*(sin(x)^2)))*(sin(x)^6))+(-1*cos(x)*(sin(x)^6))+((cos(x)+(-1*cos(x)*(sin(x)^2)))*(sin(x)^6))+(sin(x)^6))*(sin(x)^-8)))*((1+cos(x))^-1))
>>> fraction
0
>>> cos(x)/(1+sin(x))+(1+sin(x))/cos(x)-2/cos(x)
((-2*(cos(x)^-1))+((1+sin(x))*(cos(x)^-1))+(cos(x)*((1+sin(x))^-1)))
>>> fraction
0
>>> (1+1/cos(x))*cos(x)-sin(x)^2/(1-cos(x))
((-1*((1+(-1*cos(x)))^-1)*(sin(x)^2))+((1+(cos(x)^-1))*cos(x)))
>>> fraction
((1+(-1*(1+(-1*(sin(x)^2))))+(-1*(sin(x)^2)))*((1+(-1*cos(x)))^-1))
>>> fraction
0
>>> sin(x)-2*sin(x)^3-tan(x)*(2*cos(x)^3-cos(x))
((-1*((-1*cos(x))+(2*(cos(x)^3)))*(cos(x)^-1)*sin(x))+(-2*(sin(x)^3))+sin(x))
>>> fraction
(((-2*(cos(x)+(-1*cos(x)*(sin(x)^2)))*sin(x))+(-2*cos(x)*(sin(x)^3))+(2*cos(x)*sin(x)))*(cos(x)^-1))
>>> expand
0
>>> (sin(5*x)-sin(3*x))/x
(((-1*sin((3*x)))+sin((5*x)))*(x^-1))
>>> limit 0
limit x->0 (((-1*sin((3*x)))+sin((5*x)))*(x^-1)) = 2
>>> (1-cos(2*x))/x^2
((1+(-1*cos((2*x))))*(x^-2))
>>> limit 0
limit x->0 ((1+(-1*cos((2*x))))*(x^-2)) = 2
>>> 90*x^2*sin(2+6*x^3)
(90*(x^2)*sin((2+(6*(x^3)))))
>>> Sdx
(90*(((-18^-1)*cos(2)*cos((6*(x^3))))+((18^-1)*sin(2)*sin((6*(x^3))))))
>>> expand
((-5*cos(2)*cos((6*(x^3))))+(5*sin(2)*sin((6*(x^3)))))
>>> (1-cos(100*x))/sin(x)^2
((1+(-1*cos((100*x))))*(sin(x)^-2))
>>> limit 0
limit x->0 ((1+(-1*cos((100*x))))*(sin(x)^-2)) = 5000
>>> 
```
