import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

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
#conmuit()

G0 = sym.Matrix([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]])
G1 = sym.Matrix([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]])
G2 = sym.Matrix([[0,0,0,-i],[0,0,i,0],[0,i,0,0],[-i,0,0,0]])
G3 = sym.Matrix([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]])

AG = [G0,G1,G2,G3]

def anticonmutador(A,B):
    return A*B + B*A

def anticonit():
    for i in range(4):
        for j in range(4):
            if j != i:
                print((i,j), anticonmutador(AG[i],AG[j]))

ed = [1,-1,-1,-1]

M = sym.diag(*ed)

I4 = sym.eye(4)

print(2*M*I4)

#anticonit() 