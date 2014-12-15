#! /usr/bin/env python3

import numpy as np, math
gamma = np.logspace(0,4,100)
def Maxwell(gamma,avggamma):
    y = []
    for i in gamma:
        y.append( 2.0/math.sqrt(math.pi)*math.sqrt(27.0*i/(8.0*pow(avggamma,3)))*math.exp(-3.0*i/(2.0*avggamma)) )
    return y

def BrokenPL(gamma,p1,p2,gammabreak):
    y = []
    a = 1.0e5
    for i in gamma:
        if i <= gammabreak:
            y.append( a*math.pow(gammabreak,p1-p2)*pow(i,-p1) )
        elif i > gammabreak:
            y.append( a*pow(i,-p2) )
    return y

import matplotlib.pyplot as plt
plt.loglog(gamma,gamma*Maxwell(gamma,1000))
plt.loglog(gamma,gamma*BrokenPL(gamma,0.5,2.5,1000))
plt.show()

#d = Maxwell(gamma,5000)
#for i in range(len(gamma)):
#    print(gamma[i], d[i])
