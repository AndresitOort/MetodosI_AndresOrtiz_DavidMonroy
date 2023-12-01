import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import itertools as it

def GetExperiments(N=1000):
    
    freq = np.zeros(N)
    print(freq)
    for i in range(int(N)):
        
        d1 = np.random.randint(1,7)
        d2 = np.random.randint(1,7)
        
        freq[i] = d1+d2
        
    return freq

freq = GetExperiments()

print(freq)

