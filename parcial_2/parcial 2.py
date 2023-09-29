import numpy as np 
import sympy as sym
import math
import matplotlib.pyplot as plt

_x = sym.symbols('x')
N = 20

def functionnum(x):
    return (x**3)/(np.exp(x)-1)

def functionden(x):
    return (x**3)/(1-np.exp(-x))

def legendre(f,N,a,b):
    
    r,w = np.polynomial.legendre.leggauss(N)
    
    I = 0
    
    for i in range(len(r)):
         I += w[i]*(f((r[i]/2)*(b-a)+((b+a)/2)))
    
    return I*((b-a)/2)

def lagerreden(f,N):
    
    r,w = np.polynomial.laguerre.laggauss(N)
    
    I = 0
    
    for i in range(len(r)):
        
        I += w[i]*f(r[i])
    
    return I

#limites de integracion 

a = ((6.626e-34)*(3e15))/((1.3806e-23)*5772)
b = ((6.626e-34)*(7.5e14))/((1.3806e-23)*5772)

print('limites de integracion')
print(a,b)

den = lagerreden(functionden,N)

num = legendre(functionnum,N,b,a)

ans = (num/den)*100

print('porcentaje de ingreso')
print(ans)

#La diferencia entre los porcentajes se puede deber a que estamos suponiendo
#que el sol es un cuerpo negro perfecto. Sumado a esto la atmosfera tambien
#interfiere en el ingreso de rayo UV a la tierra.


def GetLaguerre(n,x):
    if n==0:
        poly = sym.Number(1)
    elif n==1:
        poly = - x + 1
    else:
        poly = ((2*(n-1)+1-x)*GetLaguerre(n-1,x)-((n-1)*GetLaguerre(n-2,x)))/(n)
   
    return sym.expand(poly,x)

print('polinomio n° 2')
print(GetLaguerre(2,_x))

def DiffLaguerre(n,x):
    pn = GetLaguerre(n,x)
    return sym.diff(pn,x,1)

def GetNewton(f,df,xn,itmax=10000,precision=1e-14):
    
    error = 1.
    it = 0
    
    while error >= precision and it < itmax:
        
        try:
            
            xn1 = xn - f(xn)/df(xn)
            
            error = np.abs(f(xn)/df(xn))
            
        except ZeroDivisionError:
            print('Zero Division')
            
        xn = xn1
        it += 1
        
    if it == itmax:
        return False
    else:
        return xn
    
def GetRoots(f,df,x,tolerancia = 10):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewton(f,df,i)

        if  type(root)!=bool:
            croot = np.round( root, tolerancia )
            
            if croot not in Roots:
                Roots = np.append(Roots, croot)
                
    Roots.sort()
    
    return Roots

def GetAllRootsGLag(n):

    j = n+((n-1)*np.sqrt(n))
    xn = np.linspace(0,j,round(50+((j**2)/2)))
    
    Legendre = np.array([])
    DLegendre = np.array([])
    
    for i in range(n+1):
        Legendre =  np.append(Legendre,GetLaguerre(i,_x))
        DLegendre = np.append(DLegendre,DiffLaguerre(i,_x))
    
    poly = sym.lambdify([_x],Legendre[n],'numpy')
    Dpoly = sym.lambdify([_x],DLegendre[n],'numpy')
    Roots = GetRoots(poly,Dpoly,xn)

    if len(Roots) != n:
        ValueError('El número de raíces debe ser igual al n del polinomio.')
    
    return Roots

def GetWeigthsLag(n,x):
    
    r = GetAllRootsGLag(n)
    
    w1 = sym.exp(-x)*((x-r[1])/(r[0]-r[1]))
    w2 = sym.exp(-x)*((x-r[0])/(r[1]-r[0]))
    
    w1i = sym.integrate(w1,(x,0,sym.oo))
    w2i = sym.integrate(w2,(x,0,sym.oo))
    
    return (w1i,w2i)
    

print('raices del 2do polinomio')
print(GetAllRootsGLag(2))

print('pesos del 2do polinomio')
print(GetWeigthsLag(2,_x))

def poli3(x):
    return(x**3)

def lagsym(f):
    r = GetAllRootsGLag(2)
    w = GetWeigthsLag(2,_x)
    
    I = 0
    
    for i in range(len(r)):
        I += w[i]*f(r[i])
        
    return I


print('resultado de integral polinomio de 3er grado')
print(lagsym(poli3))