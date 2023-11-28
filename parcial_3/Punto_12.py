import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

F = np.array([lambda x,y,z: 6*x-2*np.cos(y*z)-1,
            lambda x,y,z: 9*y+np.sqrt(x**2+np.sin(z)+1.06)+0.9,
            lambda x,y,z: 60*z+3*np.exp(-x*y)+10*np.pi-3])

G = np.array([lambda x,y: np.log(x**2+y**2)-np.sin(x*y)-np.log(2*np.pi),
              lambda x,y: np.exp(x-y)+np.cos(x*y)])


def GetF(F,r):
    
    n = r.shape[0]
    
    v = np.zeros_like(r)

    if len(F) == 3:

        for i in range(n):
            v[i] = F[i](r[0],r[1],r[2])
    else:   

        for i in range(n):
            v[i] = F[i](r[0],r[1])
            
    return v

def J(F,r,h=1e-6):

    n = r.shape[0]
    J= np.zeros((n,n))
    
    if len(F) == 3:

        for i in range(n):
            for j in range(n):

                rf= r.copy()
                rb= r.copy()

                rf[j] = rf[j]+h
                rb[j] = rb[j]-h

                J[i,j] = (F[i](rf[0],rf[1],rf[2]) - F[i](rb[0],rb[1],rb[2]))/(2*h)
    
    else: 

        for i in range(n):
            for j in range(n):

                rf= r.copy()
                rb= r.copy()

                rf[j] = rf[j]+h
                rb[j] = rb[j]-h

                J[i,j] = (F[i](rf[0],rf[1]) - F[i](rb[0],rb[1]))/(2*h)

    return J

def NewtonRaphson(F,r,itmax=100,precision=1e-7):
    
    error = 1
    it = 0
    
    while error > precision and it < itmax:
        
        rc= r 

        Fn = GetF(F,rc)
        Jn= J(F,rc)
        Jinv = np.linalg.inv(Jn)

        r1 = rc - np.dot(Jinv,Fn)
        
        error = np.max( np.abs(r1-rc) )
        
        rc = r1
        it +=1
        
    return rc

r0= np.array([0.,0.,0.])
r1= np.array([2.,2.])

A= NewtonRaphson(F,r0)
B= NewtonRaphson(G,r1)
print(A,B)