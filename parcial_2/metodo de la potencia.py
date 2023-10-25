import numpy as np
import sympy as sym

A = np.array([[1.,2.,-1.],[1.,0.,1.],[4.,-4.,5]])
b = np.array([1.,1.,1.])

def norma(v):
    return np.sqrt(np.dot(v,v.T))

def ptm(w,A):
    wt = w.T
    p = np.dot(wt,A)
    f = np.dot(p,w)
    return f

def eigenvectors(A,b,k):
    
    z = b
       
    for i in range(k):
        w = z/norma(z)  
        mu =  ptm(w,A)
        z = np.matmul(A,w)
    
    return w,mu

print(eigenvectors(A,b,10))

v,e = eigenvectors(A,b,10)

r = np.dot(A,v)
print(r)

