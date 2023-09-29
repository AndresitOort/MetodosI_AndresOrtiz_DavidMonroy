import numpy as np
import matplotlib.pyplot as plt

def funcion(x,a=0.01,R=0.5):
    return (np.sqrt((a**2)-(x**2)))/(R+x)

class inductancia:
    def __init__(self,x,f):
        
        self.x=x
        self.f=f
        self.a=0.01
        self.R=0.5
        self.y=self.f(x)
        self.h=self.x[1]-self.x[0]

        self.valorReal= np.pi*(self.R-np.sqrt((self.R**2)-(self.a**2)))
        self.integral= 0
    
    def inducTrapecio(self):

        self.integral= self.h*0.5*(self.y[0]+self.y[-1])
        suma=0

        for i in range(len(self.y[1:-1])):
            suma+=self.y[i]
        
        self.integral=self.integral+self.h*suma
        error= abs((self.integral-self.valorReal)/self.valorReal)*100
        resultado= "El valor de la inductancia fue de {0}, con un error del {1} al valor esperado".format(self.integral,error)

        return resultado

    def inducSimpson(self):
        
        h=(self.x[-1]-self.x[0])/2
        y_m= self.f((self.x[-1]+self.x[0])/2)
        self.integral=(h/3)*(self.y[1]+4*y_m+self.y[-2])
        
        error= abs((self.integral-self.valorReal)/self.valorReal)*100
        resultado= "El valor de la inductancia fue de {0}, con un error del {1} al valor esperado".format(self.integral,error)

        return resultado


x= np.linspace(-0.01,0.01,100)

induc= inductancia(x,funcion)
print(induc.inducSimpson())

y= funcion(x)
#print(x)
#plt.plot(x,y)
#plt.show()