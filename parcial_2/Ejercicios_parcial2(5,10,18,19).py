import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import math 

#Ejercicio numero 18

_x = sym.symbols('x',real=True)

def GetHermite(n,x):
    if n == 0:
        poly = sym.Number(1)
    elif n == 1:
        poly = 2*x
    else:
        poly = 2*x*GetHermite(n-1,x) - sym.diff(GetHermite(n-1,x),x,1)
    return sym.expand(poly,x)

def GetHermiteWR(n):
    
    R,W = np.polynomial.hermite.hermgauss(n)
    
    return (R,W)

def NpolinomiosHermite(n):
    
    for i in range(1,n+1):
        
        wr = GetHermiteWR(i)
        
        print('Polinomio n° ' + str(i))
        print(GetHermite(i,_x))
        print(wr[0])
        print(wr[1])
        
nh = 20

#NpolinomiosHermite(nh)

#La funcion se define por medio de la sustitucion hecha en el Imagen_1

def function(x):
    return x**4

def IntegralHermite(f,n):
    
    res = 0
    
    r,w = GetHermiteWR(n)
    
    for i in range(len(r)):
        
        res += w[i]*f(r[i])
    
    return res*(2/np.sqrt(np.pi))

print(IntegralHermite(function,20))

#Ejercicio numero 5

def GetLaguerre(n,x):
    if n==0:
        poly = sym.Number(1)
    elif n==1:
        poly = - x + 1
    else:
        poly = ((2*(n-1)+1-x)*GetLaguerre(n-1,x)-((n-1)*GetLaguerre(n-2,x)))/(n)
   
    return sym.expand(poly,x)

def DiffLaguerre(n,x):
    pn = GetLaguerre(n,x)
    return sym.diff(pn,x,1)

def GetNewton(f,df,xn,itmax=10000,precision=1e-9):
    
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
    
    weigths = np.array([])
    roots = GetAllRootsGLag(n)
    L_n1 = sym.lambdify([x],GetLaguerre(n+1,x),'numpy')
    
    for i in roots:
        
        w = i/(((n+1)**2)*(L_n1(i))**2)
        
        if w not in weigths:
            
            weigths = np.append(weigths,w)
            
    return weigths

def NpolyLaguerre(n):
    
    for i in range(1,n+1):
        print('polinomio n° ' + str(i))
        print(GetLaguerre(i,_x))
        print(GetAllRootsGLag(i))
        print(GetWeigthsLag(i,_x))

ng = 20

#NpolyLaguerre(ng)

#Ejercicio numero 19

T = np.linspace(0,20,10000)

def function19(x,T,dT):
    return (np.tanh((np.sqrt(x**2+dT**2)))*(300/(2*T)))/np.sqrt(x**2+dT**2)

"""def integral(f,T):
    
    r,w =np.polynomial.legendre.leggauss(50)
    
    for i in T:
        
    

integral(function19,25,120)    """
