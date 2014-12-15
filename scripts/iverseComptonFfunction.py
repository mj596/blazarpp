import numpy as np
import math

#define mec2  (8.18726233e-7)

def f( g, e0, e, miu ):
    w = e/g
    b = 2.0*(1.0-miu)*e0*g
    if e>e0 and e<b*g/(1+b):
        wp = 1.0-w
        fx = (1.0+0.5*w*w/wp-2.0*w/(b*wp)+2.0*w*w/(b*b*wp*wp))
        if fx>0.0:
            return fx
        else:
            return 0.0
    else:
        return 0.0

def beta( _g ):
    return math.sqrt( math.pow(_g,2) - 1.0 )/_g

eV2erg = 1.60217646e-12
mec2 = 8.18726233e-7
PLANCK_H = 6.626e-27
SIGMA_T = 6.652453e-25
vext = 10.0*eV2erg/PLANCK_H


gamma = np.logspace(2,2,1)
v = np.logspace(15,30,1000)
theta = 0.0667
Gamma = 15.0

e0 = vext*PLANCK_H*Gamma/mec2
miu = -( math.cos(theta)-beta(Gamma) )/( 1.0-beta(Gamma)*math.cos(theta) )


import matplotlib.pyplot as plt
for g in gamma:
    print g
    fvec = []
    for i in v:
        e = PLANCK_H*i/mec2
        fvec.append( f(g,e0,e,miu) )
    plt.loglog( v, fvec )
plt.show()
