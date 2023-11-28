import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

# Punto 4

def probBirthday(n):

    prob=1
    problist=np.array([])
    N=np.linspace(0,n,n)

    for i in range(n):
        
        prob*=(365-i)/(365)
        problist=np.append(problist,prob)

    plt.scatter(N,problist)
    plt.show()

    return #problist[-1]

print(probBirthday(80))

#Punto 8

