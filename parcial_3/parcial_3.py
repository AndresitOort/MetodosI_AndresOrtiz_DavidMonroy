import numpy as np
import sympy as sym

A = np.array([[0.2,0.1,1.,1.,0.],[0.1,4.,-1.,1.,-1.],[1.,-1.,60.,0.,-2.],[1.,1.,0.,8.,4.],[0.,-1.,-2.,4.,700.]])
b = np.array([1.,2.,3.,4.,5.])
x0 = np.array([1.,1.,1.,1.,1.])
def productointerno(Vt,A,v):
    p = np.dot(Vt,A)
    return np.dot(p,v)

def DecensoConjugado(A,b,x0,epsilon=0.01):
    
    r0 = np.dot(A,x0)-b
    p0 = -r0
    k = 0
    
    while np.max(r0) > epsilon:
        
        alf = -(np.dot(r0.T,p0)/productointerno(p0.T,A,p0))
        x0 = x0 + alf*p0
        r0 = np.dot(A,x0) - b
        bet = productointerno(r0.T,A,p0)/productointerno(p0.T,A,p0)
        p0 = -r0 + bet*p0
        
        k += 1
    
    return x0,k

print(DecensoConjugado(A,b,x0))