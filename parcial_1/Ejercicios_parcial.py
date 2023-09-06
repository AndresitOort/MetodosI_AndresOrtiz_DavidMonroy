import numpy as np
import sympy as sp

x = np.linspace(0.56,58,100)

l = np.zeros([3,3])

x_0 = 0.5

x_1 = 0.6

print(x_0*x_1)

def function(x):
    return (np.e)**(-x)-x

def derivate(f,x,h=1e-6):
    return(f(x+h)-f(x-h))/(2*h)

def newton_gregory(f, df):
    return None

#def newton_greogory(x):
    