import numpy as np 
import sympy as sym
import math
import matplotlib.pyplot as plt

_x = sym.symbols('x')
N = 20

def functionnum(x):
    return (x**3)/(np.exp(x)-1)

def functionden(x):
    return (x**3)/(1-np.exp(-x))

def legendre(f,N,a,b):
    
    r,w = np.polynomial.legendre.leggauss(N)
    
    I = 0
    
    for i in range(len(r)):
         I += w[i]*(f((r[i]/2)*(b-a)+((b+a)/2)))
    
    return I*((b-a)/2)

def lagerreden(f,N):
    
    r,w = np.polynomial.laguerre.laggauss(N)
    
    I = 0
    
    for i in range(len(r)):
        
        I += w[i]*f(r[i])
    
    return I

#limites de integracion 

a = ((6.626e-34)*(3e15))/((1.3806e-23)*5772)
b = ((6.626e-34)*(7.5e14))/((1.3806e-23)*5772)
print(a,b)

den = lagerreden(functionden,N)

num = legendre(functionnum,N,b,a)

ans = (num/den)*100

print(ans)

#La diferencia entre los porcentajes se puede deber a que estamos suponiendo
#que el sol es un cuerpo negro perfecto.

