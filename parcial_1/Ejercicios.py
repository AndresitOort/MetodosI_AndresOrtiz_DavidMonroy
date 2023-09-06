import numpy as np
import matplotlib.pyplot as plt
import math as m
import sympy as sp

x = np.linspace(0.1,1.1,100)

def function(x):
    return np.sqrt(np.tan(x))

def derivada_p(f,x,h=0.01):
    return(f(x+h)-f(x))/h

def derivada_c(f,x,h=0.01):
    return(f(x+h)-f(x-h))/(2*h)

def derivada_e(x):
    return ((1/np.cos(x))**2)/(2*np.sqrt(np.tan(x)))

y_0 = function(x)
y_1 = derivada_p(function,x)
y_2 = derivada_c(function,x)
y_3 = derivada_e(x)

def error_nodal(der,f,x):
    return abs(der(f,x)-derivada_e(x))

err_p = error_nodal(derivada_p,function,x)
err_c = error_nodal(derivada_c,function,x)

plt.scatter(x,err_p)
plt.scatter(x,err_c)
plt.show()

plt.scatter(x,y_0,label='funcion')
plt.scatter(x,y_1,label='derivada progresiva')
plt.scatter(x,y_2,label='derivada central')
plt.plot(x,y_3,label='derivada exacta')
plt.legend()
plt.show()