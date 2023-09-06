import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sym
import math

x = np.linspace(0,6,100)

def Lagrange(x,X,i):
    L = 1
    for j in range(X.shape[0]):
        if i != j:
            L *= (x - X[j])/(X[i]-X[j])
    return L

def Interpolate(x,X,Y):
    Poly = 0
    for i in range(X.shape[0]):
        Poly += Lagrange(x,X,i)*Y[i]
    return Poly

def archivo_leer(ruta):
    g = -9.8
    _x = sym.Symbol('x',real=True)
    a = open(ruta)
    ejes = a.readline().split(",")
    x_1=np.array([])
    y_1=np.array([])
    L = a.readline()
    while L != '':
        l = L.split(",")
        x_1 = np.append(x_1,float(l[0]))
        y_1 = np.append(y_1,float(l[1]))
        L = a.readline()   
    y_2 = Interpolate(_x,x_1,y_1)
    _y_2= sym.simplify(y_2)
    n = sym.poly_from_expr(_y_2)
    coeffs = n[0].all_coeffs()
    angle_r = math.atan(coeffs[1])
    angle_deg = round(angle_r*180/np.pi)
    print('el angulo de salida es: ' + str(angle_deg)+ '°')
    v_0 = round(math.sqrt(g/(2*coeffs[0]*(np.cos(angle_r)**2))))
    print('la velocidad inicial es: ' + str(v_0) + 'm/s^2')

archivo_leer("C:\\Users\\david\\OneDrive\\Documentos\\programas\\Metodos\\MetodosComputacionales202320\\Parabolico.csv")


