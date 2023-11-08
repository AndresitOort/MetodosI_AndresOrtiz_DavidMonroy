import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

#PUNTO NUMERO 4

print("PUNTO NUMERO 4")

A = np.array([[1.,0.,0.],[5.,1.,0.],[-2.,3.,1.]])
B = np.array([[4.,-2.,1.],[0.,3.,7.],[0.,0.,2.]])

def MultiMatAux(A,B):
    
    n,m = (A.shape[0],B.shape[1])
        
    D = np.zeros((n,m))
    
    for i in range(n):
        for j in range(m):
            l = sum(A[i,:]*B[:,j])
            D[i,j] = l
            
    return D

def MultiMat(A,B):     
    if len(B.shape) == 1:
        
        V = B.reshape(B.shape[0],1)
        
        if A.shape[1] != V.shape[0]:
            return 'No cumplen el requisito de dimensiones'
        
        D = MultiMatAux(A,V)
        
    else:
        
        if A.shape[1] != B.shape[0]:
            return 'No cumplen el requisito de dimensiones'
        
        D = MultiMatAux(A,B)
    
    return D

print(MultiMat(A,B))

#PUNTO NUMERO 7

print("PUNTO NUMERO 7")

def SOR(A,b,x0,w,tol=1e-6,itmax=1000):
    
    x = x0.copy()
    u = x.copy()
    sumk1 = x.copy()
    sumk2 = x.copy()
    
    residuo = np.linalg.norm( np.dot(A,x) - b)
    
    it = 0
    
    while it < itmax and residuo > tol:
        
        u[:] = 0
        sumk1[:] = 0
        sumk2[:] = 0
        
        for i in range(A.shape[0]):
            if i >= 1:
                for j in range(0,i):
                    sumk1[i] += A[i,j]*u[j]
                
            for j2 in range(i+1,A.shape[0]):
                sumk2[i] += A[i,j2]*x[j2]
            
            u[i] = (1-w)*x[i] + (w/A[i,i])*(b[i]-sumk1[i]-sumk2[i])
        
        x = u.copy()
        
        residuo = np.max(np.abs(np.dot(A,x) - b))
        
        it += 1
    
    return x,it

def FindMinItSOR(A,b,x0,N=100):
    
    x = np.linspace(0.5,1.5,N)
    ini = N
    
    omega = 0
    
    for i in x[1:-1]:
        
        v,it = SOR(A,b,x0,i)
        
        if it < ini:
            ini = it
            omega = i
    
    return SOR(A,b,x0,omega),omega

A1 = np.array([[3,-1,-1],[-1.,3.,1.],[2,1,4]])
b1 = np.array([1.,3.,7.])

x0 = np.array([0.,0.,0.])

print(FindMinItSOR(A1,b1,x0))

#PUNTO NUMERO 15

print("PUNTO NUMERO 15")

i = sym.I

MAx = sym.Matrix([[0,1],[1,0]])
MAy = sym.Matrix([[0,-i],[i,0]])
MAz = sym.Matrix([[1,0],[0,-1]])

AM = [MAx,MAy,MAz]

def conmutador(A,B):
    return A*B - B*A

def conmuit():
    for i in range(3):
        for j in range(3):
            if j != i:
                print((i,j), conmutador(AM[i],AM[j]))

conmuit()

#PUNTO NUMERO 16

print("PUNTO NUMERO 16")

G0 = sym.Matrix([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]])
G1 = sym.Matrix([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]])
G2 = sym.Matrix([[0,0,0,-i],[0,0,i,0],[0,i,0,0],[-i,0,0,0]])
G3 = sym.Matrix([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]])

AG = [G0,G1,G2,G3]

def anticonmutador(A,B):
    return A*B - B*A

def anticonit():
    for i in range(4):
        for j in range(4):
            if j != i:
                print((i,j), anticonmutador(AG[i],AG[j]))

anticonit() 