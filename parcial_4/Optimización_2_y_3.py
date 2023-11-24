import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import sympy as sym

#Punto 2

x,y,z= sym.symbols('x,y,z',real=True)

    # Definimos la función a optimizar

def f(p):
    return p[0]**2+p[1]**2+p[2]**2-2*p[2]+1

    # Definimos las restricciones bajo las que se optimiza y elegimos un punto semilla

restriccion= ( {'type':'eq','fun': lambda p: 2*p[0]-4*p[1]+5*p[2]-2})
p0= [3,4,1]

minimo= spo.minimize(f,p0,constraints=restriccion)

#print("El mínimo de la función bajo la restricción es:{0} ".format(minimo.fun))

#Punto 3

 # Definimos la función de volumen

def V(p):
    return -(p[0]*p[1]*p[2])

    # Definimos las restricciones para el volumen de la caja

area= ( {'type':'eq','fun': lambda p: p[0]*p[1]+2*p[1]*p[2]+2*p[0]*p[2]-12})
pi=[2,3,4]

volumen = spo.minimize(V,pi,constraints=area)

print("El máximo volumen para una caja con área lateral de 12cm es: {0} cm^3".format(round(abs(volumen.fun))))
