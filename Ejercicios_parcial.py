import numpy as np
import sympy as sp

x = np.linspace(0.56,58,100)

def boltzano(f, x0, x1, tol=1000):
    
    x_0 = x0
    x_1 = x1
    
    i = 0
    
    while tol >= i and f(x_0)*f(x_1) >= 1e-10:
        
        x_0 = x_0/2
        x_1 = x_1/2
        
        i += 1
    
    return(x_0,x_1)

def pm(f, t):
    
    x2 = (t[1]+t[0])/2
    
    f2=f(x2)
    
    return np.array([t[0],x2,t[1]])

def function(x):
    return (np.e)**(-x)-x

t_x = boltzano(function,0,1)
a=pm(function,t_x)

def derivate(f,x,h=1e-6):
    return(f(x+h)-f(x-h))/(2*h)

def muller(f,a, tol=100):
    
    l = np.zeros([3,3])
    
    l[:,0]=f(a)
    
    d=np.array([])
    
    for i in range(1,len(a)):
        for j in range(i,len(a)):
            l[j,i] = l[j,i-1]-l[j-1,i-1]
            
            
    d = np.append(d,l[0,0])
    d = np.append(d,l[1,1])
    d = np.append(d,l[2,2])
    
    return d

d=muller(function,a)

def buscar_x3(d,a):
    a_ = d[2]
    b_ = d[1]-(a[0]+a[1])*a_
    c_ = d[0]-(a[0]*d[1])+(a[0]*a[1]*a_)
    
    if b_ < 0:
        x_3 = (-2*c_)/(b_-np.sqrt(b_**2-4*a_*c_))
    else:
        x_3 = (-2*c_)/(b_+np.sqrt(b_**2-4*a_*c_))
    
    return x_3

x_3 = buscar_x3(d,a)

def Newton(f,df,xn,itmax=100):
    
    precision=1e-8
    
    error=1.
    it=0

    while error>precision and it<itmax:

        try:
            xn1= xn - f(xn)/df(f,xn)
            error = np.abs(f(xn)/df(f,xn))
       
        except ZeroDivisionError:
            print("División por cero :(")
       
        xn=xn1
        it+=1

    if it == itmax:
        return False
    else:
        return xn

print(Newton(function,derivate,x_3))

#error>precision
print(t_x,a)