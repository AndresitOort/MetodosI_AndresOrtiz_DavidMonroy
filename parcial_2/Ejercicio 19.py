#Ejercicio numero 19

T = np.arange(1,20,10e-4)

def function19(x,T,DT):
    return (sym.tanh((sym.sqrt(x**2+dT**2)))*(300/(2*T)))/sym.sqrt(x**2+dT**2)

def integral(f,T):
    
    r,w =np.polynomial.legendre.leggauss(50)
    
    for i in T:
        
        
        
integral(function19,T)
