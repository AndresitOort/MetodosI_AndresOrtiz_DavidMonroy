import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

x = sym.Symbol('x',real=True)
y = sym.Symbol('y',real=True)

def f(x,y):
    
    z = x + sym.I*y
    f = z**3 - 1
    f = f.expand()

    return [sym.re(f),sym.im(f)]

#print(f(x,y))

def J(f):

    J = sym.zeros(2,2)

    for i in range(2):
        for j in range(2):
            if j==0:
                J[i,j] = sym.diff(f[i],x,1)
            else:
                J[i,j] = sym.diff(f[i],y,1)
    
    return J

def NewtonRaphsonI(z,Fn,Jn,itmax=100,precision=1e-7):
    
    error = 1
    it = 0
    
    while error > precision and it < itmax:
        
        IFn = Fn(z[0],z[1])
        IJn = Jn(z[0],z[1])
        
        z1 = z - np.dot(IJn,IFn)
        
        error = np.max( np.abs(z1-z) )
        
        z = z1
        it +=1
        
    return z

def GetFractal(Fn,Jn,N=500):

    x=np.linspace(-1,1,N)
    y=np.linspace(-1,1,N)
    Fractal = np.zeros((N,N), np.int64)
    r=0

    for i in range(len(x)):
        for j in range(len(y)):
            
            z= np.array([x[i],y[j]])

            r= NewtonRaphsonI(z,Fn,Jn)

            if r[0]>0:
                Fractal[i][j]= 255 
            elif r[1]<0:
                Fractal[i][j]= 100
            else:
                Fractal[i][j]= 20
            
    
    plt.imshow(Fractal, cmap='coolwarm' ,extent=[-1,1,-1,1])
    plt.show()
    
    return "Bonito, ¿no?"

A= f(x,y)
B= J(A).inv()

An= sym.lambdify ([x,y],A,'numpy') 
Bn= sym.lambdify ([x,y],B,'numpy')

print(GetFractal(An,Bn))