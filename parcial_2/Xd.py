import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import math

_x=sym.symbols('x')

def funcion(x):
    return x**2

def Legendre(x,n):

    polys=np.array

    for i in range(n):
        base= (x**2 - 1)**i
        diff= sym.diff(base,x,i)
        fact= 1/((2**i) * (sym.factorial(i)))
        print (diff*fact)


print(Legendre(_x,5))