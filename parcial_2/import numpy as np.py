import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import math

A = np.array([[2,1,1],\
     [1,1,-2],
     [5,10,5]])

A1 = np.array([[1,2,1,-1],\
             [3,2,4,4],
             [4,4,3,4],
             [2,0,1,5]])

b = np.array([8,-2,10])

n = np.shape(A)[0]
M = np.zeros(shape=(n,n+1))

M[:,0:n] = A
M[:,n] = b

def back_substitution(A, b):
    '''
    Args:
        A (np.array): Matriz de coeficientes triangular superior
        b (np.arra): vector de constantes
        
    '''

    n = np.shape(A)[0]

    x = np.zeros(n)
    for i in range(n-1,-1,-1):
        sum = b[i]
        for j in range(n-1,i,-1):
            print(i,j)
            sum -= A[i,j]*x[j]
        x[i] = sum/A[i,i]

    return x

def leer_m(A):
    
    n = np.shape(A)[0]
    
    c = np.array([])
    
    for i in range(1,n):
        for j in range(i):
            v = A[i,j]/A[j,j]
            c = np.append(c,v)
            
    return c

def gaussian_elimination(A,b):
    
    n = np.shape(A)[0]
    
    D = np.zeros(shape=(n,n+1))
    
    D[:,0:n] = A
    D[:,n] = b
    
    C = leer_m(A)
    t = 0
    for j in range(n-1):
        for i in range(j+1,n):
            D[i,:] = D[i,:] - D[j,:]*C[t]
            C = leer_m(D)
            t += 1
        
    return D

print(leer_m(A))
print(gaussian_elimination(A,b))