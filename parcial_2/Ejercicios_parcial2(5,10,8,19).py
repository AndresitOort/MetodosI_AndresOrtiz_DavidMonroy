import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import math 

_x = sym.symbols('x')

def GetNewton(f,df,xn,itmax=10000,precision=1e-12):
    
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
    
def GetRoots(f,df,x,tolerancia = 5):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewton(f,df,i)

        if  type(root)!=bool:
            croot = np.round( root, tolerancia )
            
            if croot not in Roots:
                Roots = np.append(Roots, croot)
                
    Roots.sort()
    
    return Roots

def GetHermite(n,x):
    if n == 0:
        poly = sym.Number(1)
    elif n == 1:
        poly = 2*x
    else:
        poly = 2*x*GetHermite(n-1,x) - sym.diff(GetHermite(n-1,x),x,1)
    return sym.expand(poly,x)

def DiffHermite(n,x):
    pn = GetHermite(n,x)
    return sym.diff(pn,x,1)

def GetAllRootsGHer(n):
    
    c = np.sqrt(4*n+1)
    xn = np.linspace(-c,c,100)
    
    Hermite = np.array([])
    DHermite = np.array([])
    
    for i in range(n+1):
        Hermite =  np.append(Hermite,GetHermite(i,_x))
        DHermite = np.append(DHermite,DiffHermite(i,_x))
    
    poly = sym.lambdify([_x],Hermite[n],'numpy')
    Dpoly = sym.lambdify([_x],DHermite[n],'numpy')
    Roots = GetRoots(poly,Dpoly,xn)

    if len(Roots) != n:
        ValueError('El número de raíces debe ser igual al n del polinomio.')
    
    return Roots

def GetWeigthsHer(n,x):
    
    roots = GetAllRootsGHer(n)
    l_h1 = sym.lambdify([x],GetHermite(n-1,x),'numpy')
    weigths = np.array([])
    
    for i in roots:
        
        w = (2**(n-1)*math.factorial(n)*np.sqrt(np.pi))/((n**2)*(l_h1(i)**2))
            
        weigths = np.append(weigths,w)
    
    return weigths

def NpolinomiosHermite(n):
    
    for i in range(n+1):
        if n == 0:
            print('El polinomio 0 no tiene raices')
        else:
            print('polinomion n° '+str(i))
            print(GetHermite(i,_x))
            print(GetAllRootsGHer(i))
            print(GetWeigthsHer(i,_x))

nh = 20

#NpolinomiosHermite(nh)

#La funcion se define por medio de la sustitucion hecha en el Imagen_1

def function(x):
    return x**4

def IntegralHermite(f,n):
    
    res = 0
    
    r = GetAllRootsGHer(n)
    w = GetWeigthsHer(n,_x)
    
    for i in range(len(r)):
        
        res += w[i]*f(r[i])
    
    return res*(2/np.sqrt(np.pi))



print(IntegralHermite(function,20))