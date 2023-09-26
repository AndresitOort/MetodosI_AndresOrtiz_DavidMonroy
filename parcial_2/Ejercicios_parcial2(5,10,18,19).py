import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import math 

#Ejercicio numero 18

_x = sym.symbols('x')

def GetHermite(n,x):
    if n == 0:
        poly = sym.Number(1)
    elif n == 1:
        poly = 2*x
    else:
        poly = 2*x*GetHermite(n-1,x) - sym.diff(GetHermite(n-1,x),x,1)
    return sym.expand(poly,x)

def GetHermiteWR(n):
    
    R,W = np.polynomial.hermite.hermgauss(n)
    
    return (R,W)

def NpolinomiosHermite(n):
    
    for i in range(1,n+1):
        
        wr = GetHermiteWR(i)
        
        print('Polinomio n° ' + str(i))
        print(GetHermite(i,_x))
        print(wr[0])
        print(wr[1])
        
nh = 20

#NpolinomiosHermite(nh)

#La funcion se define por medio de la sustitucion hecha en el Imagen_1

def function(x):
    return x**4

def IntegralHermite(f,n):
    
    res = 0
    
    r,w = GetHermiteWR(n)
    
    for i in range(len(r)):
        
        res += w[i]*f(r[i])
    
    return res*(2/np.sqrt(np.pi))

print(IntegralHermite(function,20))

#Ejercicio numero 5

def GetLaguerre(n,x):
    if n==0:
        poly = sym.Number(1)
    elif n==1:
        poly = - x + 1
    else:
        poly = ((2*(n-1)+1-x)*GetLaguerre(n-1,x)-((n-1)*GetLaguerre(n-2,x)))/(n)
   
    return sym.expand(poly,x)

def GetLaguerreWR(n):
    
    Roots,Weihts = np.polynomial.laguerre.laggauss(n)
    
    return (Roots,Weihts)

def NpolyLaguerre(n):
    
    for i in range(1,n+1):
        wr = GetLaguerreWR(i)
        print('polinomio n° ' + str(i))
        print(GetLaguerre(i,_x))
        print(wr[0])
        print(wr[1])

ng = 20

#NpolyLaguerre(ng)

#Ejercicio numero 19

