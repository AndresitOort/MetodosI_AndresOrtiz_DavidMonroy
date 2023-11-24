import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

C= 729/10000
x = sym.Symbol('x',real=True)
y = sym.Symbol('y',real=True)
r= np.array([80.])

def f(x,y):
    
    z = sym.sin(x+sym.I*y)
    f = z**6 + C*z**2 - C
    f = f.expand()
    return sym.re(f),sym.im(f)

k1,k2 = f(x,y)

print(k1,k2)

def J(F,r,h=1e-6):

    n = 2
    J= np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):

            rf= r.copy()
            rb= r.copy()

            rf[j] = rf[j]+h
            rb[j] = rb[j]-h

            J[i,j] = (F[i](rf[0]) - F[i](rb[0]))/(2*h)

    return J

def NewtonRaphson(F,r,itmax=100,precision=1e-7):
    
    error = 1
    it = 0
    
    while error > precision and it < itmax:
        
        rc= r 

        Fn = F(F,rc)
        Jn= J(F,rc)
        Jinv = np.linalg.inv(Jn)

        r1 = rc - np.dot(Jinv,Fn)
        
        error = np.max( np.abs(r1-rc) )
        
        rc = r1
        it +=1
        
    return rc

def GetAllRoots(F,x,e=10):
    
    Roots = np.array([])
    
    for i in x:

        root=NewtonRaphson(F,i)
        
        croot = np.round(root, e)
        
        if croot not in Roots:
            Roots = np.append(Roots,croot)
            
    Roots.sort()
    
    return Roots

'''x= np.linspace(0,2*np.pi,100)

K= GetAllRoots(f,x)

print(K)'''