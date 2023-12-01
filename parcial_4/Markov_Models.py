import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from itertools import permutations
from itertools import combinations
from itertools import combinations_with_replacement

States = np.array([0,1])
#Prior =  np.array([0.2,0.8])
Prior =  np.array([0.0,0.2])

E = np.array([[0.5,0.9],[0.5,0.1]])
T = np.array([[0.8,0.2],[0.2,0.8]])

DictM = {0:'Justa',1:'Sesgada'} 
DictO = {0:'C',1:'S'}
Obs = np.array([1,0,0,0,1,0,1,0])

def GetStates(States,N=8):
    
    CStates = list( combinations_with_replacement(States,N) )
    
    #print(CStates)
    
    Permu = []
    
    for it in CStates:
        p = list(permutations(it,N))
        
        for i in p:
            if i not in Permu:
                Permu.append(i)
                
    return np.array(Permu)

HiddenStates = GetStates(States)

def GetProb(T,E,Obs,State,Prior):
    
    n = len(Obs)
    p = 1.
    
    p *= Prior[State[0]]
    
    # Matriz de transicion
    for i in range(n-1):
        p *= T[ State[i+1], State[i]  ]
        
    for i in range(n):
        p *= E[ Obs[i], State[i] ]
        
    return p

Nobs = HiddenStates.shape[0]

JObs = np.zeros(Nobs)

for j in range(Nobs):
    
    dim = HiddenStates.shape[0]
    J = np.zeros(dim)
    
    for i in range(dim):
        J[i] = GetProb(T,E,HiddenStates[j],HiddenStates[i],Prior)
        
    JObs[j] = np.sum(J)

P = np.zeros(HiddenStates.shape[0], dtype=np.float64)

for i in range(P.shape[0]):
    P[i] = GetProb(T,E,Obs,HiddenStates[i],Prior)

ii = np.where( P == np.amax(P))

print("La suma de probabilidades de todos los estados observables es: {}".format(sum(JObs)))

plt.plot(P,color='black',label='Probabilidad por secuencia')
plt.scatter(ii,max(P),color='r',label='MaxP')
plt.legend(loc=0)
plt.show()
