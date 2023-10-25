import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

A1 = np.array([[1.,0.,0.],[5.,1.,0.],[-2.,3.,1.]])
B1 = np.array([[4.,-2.,1.]])

print(A1,B1)
def MultiMat(A,B):
    
    n,m = (B.shape[0],A.shape[1])
    
    D = np.zeros((m,n))
    
    for i in range(m):
        for j in range(n):
            l = sum(A[i]*B[:,j])
            D[i,j] = l

    
    return D

print(MultiMat(A1,B1))



A = np.array([[3.,-1.,-1.],[-1.,3.,1],[2.,1.,4.]])
c = np.array([1.,3.,7.])

def norma(v):
    return np.sqrt(np.dot(v,v.T))

def SOR(A,b,w,tol=1e-6,itmax=1000):
    
    n = b.shape[0]
    
    x = np.zeros(n)
    
    err = 1
    it = 0
    
    while it < itmax and err > tol:
        
        xn = x.copy()
        
        for i in range(n):
            
            xn[i] = (1-w)*x[i]
            sum1 = 0
            sum2 = 0
            for j in range(1,i):
                sum1 += A[i,j]*x[j]
                for k in range(i+1,n):
                    sum2 += A[i,k]*x[k]
                    xn[i] += (w/A[i,i])*(b[i]-sum1-sum2)
        
        x = xn
        
        err = norma(b - np.dot(A,x))
        it += 1
    
    return x
    
print(SOR(A,c,1.5))
