import numpy as np
import math
import matplotlib.pyplot

N=1000
R=1

Volumen=np.array([])

def Volumen(R,N):

    x=np.linspace(-R,R,N)
    y=np.linspace(-R,R,N)

    area = (x[1]-x[0])*(y[1]-y[0])

    volumen = 0
    valor = 0 
    
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            
            if np.sqrt(x[i]**2+y[j]**2) > R:
                valor+=0
            else:
                valor += np.sqrt((-1*(x[i])**2)+(-1*(y[j])**2)+R**2)
            
            if np.sqrt(x[i+1]**2+y[j]**2) > R:
                valor+=0
            else:
                valor += np.sqrt((-1*(x[i+1])**2)+(-1*(y[j])**2)+R**2)

            if np.sqrt(x[i]**2+y[j+1]**2) > R:
                valor+=0
            else:
                valor += np.sqrt((-1*(x[i])**2)+(-1*(y[j+1])**2)+R**2)
             
            if np.sqrt(x[i+1]**2+y[j+1]**2) > R:
                valor+=0
            else:
                valor += np.sqrt((-1*(x[i+1])**2)+(-1*(y[j+1])**2)+R**2) 
                    
            volumen += (valor/4)*area
            valor = 0


    return volumen

print(Volumen(1,N))