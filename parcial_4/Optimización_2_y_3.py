import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import sympy as sym

#Punto 2

x,y,z= sym.symbols('x,y,z',real=True)

    # Definimos la función a optimizar

def f(p):
    return p[0]**2+p[1]**2+p[2]**2-2*p[2]+1

    # Definimos las restricciones bajo las que se optimiza y elegimos un punto semilla

restriccion= ( {'type':'eq','fun': lambda p: 2*p[0]-4*p[1]+5*p[2]-2})
p0= [3,4,1]

minimo= spo.minimize(f,p0,constraints=restriccion)

#print("El mínimo de la función bajo la restricción es:{0} ".format(minimo.fun))

#Punto 3

 # Definimos la función de volumen

def V(p):
    return -(p[0]*p[1]*p[2])

    # Definimos las restricciones para el volumen de la caja

area= ( {'type':'eq','fun': lambda p: p[0]*p[1]+2*p[1]*p[2]+2*p[0]*p[2]-12})
pi=[2,3,4]

volumen = spo.minimize(V,pi,constraints=area)

#print("El máximo volumen para una caja con área lateral de 12cm es: {0} cm^3".format(round(abs(volumen.fun))))

#Por multiplicadores de lagrange:

Lag= np.array([lambda x,y,z,p: y*z - p*y - p*2*z, 
               lambda x,y,z,p: x*z - p*x - p*2*z, 
               lambda x,y,z,p: x*y - p*2*x+p*2*y,
               lambda x,y,z,p: p*0+x*y+2*x*z+2*y*z-12])

def GetF(Lag,r):
    
    n = r.shape[0]
    
    v = np.zeros_like(r)
    
    for i in range(n):
        v[i] = Lag[i](r[0],r[1],r[2],r[3])
        
    return v

def GetJacobian(f,r,h=1e-6):
    
    n = r.shape[0]
    
    J = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            
            rf = r.copy()
            rb = r.copy()
            
            rf[j] = rf[j] + h
            rb[j] = rb[j] - h
            
            J[i,j] = ( f[i](rf[0],rf[1],rf[2],rf[3]) - f[i](rb[0],rb[1],rb[2],rb[3])  )/(2*h)
            
    
    return J

def NewtonRaphson(G,r,itmax=1000,error=1e-9):
    
    it = 0
    d = 1.
    dvector = []
    
    while d > error and it < itmax:
        
        # Vector actual
        rc = r
        
        F = GetF(G,rc)
        J = GetJacobian(G,rc)
        InvJ = np.linalg.inv(J)
        
        r = rc - np.dot(InvJ,F)
        
        diff = r - rc
        
        d = np.max( np.abs(diff) )
        
        dvector.append(d)
        print(dvector)
        
        it += 1
    
    print(it)
    return r,dvector

v = GetF(Lag,np.array([1,1,1,1]))
d= GetJacobian (Lag,np.array([1,1,1,1]))
print(d)