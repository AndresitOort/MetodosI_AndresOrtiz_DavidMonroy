import numpy as np
import sympy as sym

x = sym.symbols('x')
y = sym.symbols('y')

f = (2/3) * (x+2*y)

#Para verificar si la primera parte de la función (2/3)*(x+2y) es positiva en el dominio de x y y se realiza el siguiente despeje (2/3)*(x+2y) >= 0 -->  (x+2y) >= 0. Al evaluar en el dominio [0,1] * [0,1] la función será positiva en todo momento.
#Por lo tanto cumple la primera condicion para ser una distribuion continua de probabilidad 


total = sym.integrate(f,(x,0,1),(y,0,1))
print(total)

r1 = sym.integrate(f,(y,0,1))
print("La distribución marginal de g(x) es {}".format(r1))

r2 = sym.integrate(f,(x,0,1))
print("La distribución marginal de h(y) es {}".format(r2))

re1 = sym.integrate((x*r1),(x,0,1))
print("La media del valor esperado de X es {}".format(re1))

re2 = sym.integrate((y*r2),(y,0,1))
print("La media del valor esperado de Y es {}".format(re2))

Exy = sym.integrate((f*x*y),(x,0,1),(y,0,1))

c_1 = Exy - (re1*re2)
print("La covarianza es: {}".format(c_1))

c_2 = sym.integrate(f*(x-re1)*(y-er2),(x,0,1),(y,0,1))
print("La covarianza es: {}".format(c_2))

#Dado que la covarianza es distinta de 0 se puede determinar que las variables x y y son de pendientes
#la una de la otra
