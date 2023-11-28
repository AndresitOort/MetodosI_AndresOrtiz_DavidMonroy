import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

u1= np.array([3.,1.,0.,1.])
u2= np.array([1.,2.,1.,1.])
u3= np.array([-1.,0.,2.,-1.])
bs= np.array([-3.,-3.,8.,9.])
Vectores = np.array([u1,u2,u3])

#Base Ortonormal

def Ortonormal(Vectors):

    suma = 0
    Vectors[0] = Vectors[0]/np.linalg.norm(Vectors[0])

    for k in range(1,len(Vectors)):
    
        for j in range(k):

            suma += ((np.dot(Vectors[k],Vectors[j]))/(np.linalg.norm(Vectors[j])**2))*Vectors[j]
        
        Ort = Vectors[k]-suma
        Vectors[k] = Ort/np.linalg.norm(Ort)
        suma=0

    return Vectors

def SolOrtogonal(u1,u2,u3,bs):
    
    gen= np.array([u1,u2,u3])
    A = np.ones((4,4))
    
    for i in range(1,4):
        A[:,i] = gen[i-1]**i
        
    AT = np.dot(A.T,A)
    bT = np.dot(A.T,bs)

    xsol = np.linalg.solve(AT,bT)

    d = np.dot(A,xsol)-bs
    
    return d

#Vectores Ortonormales
O = Ortonormal(Vectores)
v0 = O[0]
v1 = O[1]
v2 = O[2]

# Con los vectores u1,u2,u3:

Sol = SolOrtogonal(u1,u2,u3,bs)
print(Sol)

# Con los vectores v0, v1, v2:

SolOrt = SolOrtogonal(v0,v1,v2,bs)
print(SolOrt)




