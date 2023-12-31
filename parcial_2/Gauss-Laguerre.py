import numpy as np
import matplotlib.pyplot as plt
import math

def funcion(x):
    return((x**3)/((np.exp(x))-1))

def GaussLaguerre(f,n):
    
    I=0
    
    error=np.array([])
    it = np.array([])

    for i in range(2,10):
        Root,Weigh=np.polynomial.laguerre.laggauss(i)
        for j in range(len(Root)):
            I+=Weigh[j]*f(Root[j])
        error=np.append(error,I/((np.pi**4)/15))
        it=np.append(it,i)
        I=0
    
    plt.scatter(it,error)
    plt.show()
    

    for i in range(n):
        R,W=np.polynomial.laguerre.laggauss(n)
        I+=W[i]*f(R[i])

    return I

print(GaussLaguerre(funcion,15))
