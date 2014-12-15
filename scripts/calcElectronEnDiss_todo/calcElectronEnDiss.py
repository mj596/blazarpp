import numpy as np
import sys
import math
import matplotlib.pyplot as plt

def read_data( filename ):
    x=[]
    y=[]
    file = open(filename,'r')
    for line in file.readlines():
        x.append(float(line.split(" ")[0]))
        y.append(float(line.split(" ")[1]))
        
        file.close()
    return np.array([x,y])

#def getEff( idin, idout ):
#
#    ein = read_data( 'Injection_'+str(id) )
#    eout = read_data( 'Ngamma_'+str(id) )
#    ein[1] *= ein[0]*ein[0]
#    eout[1] *= eout[0]*eout[0]
#
#    def calcEff( xin, xout ):
#        return np.trapz( xout[1] )/np.trapz( xin[1] )
##        return np.trapz( xout[1], x=xout[0] )/np.trapz( xin[1], x=xin[0] )
#
#    return calcEff( ein, eout )
#
#def plot_eff( ):
#    import os
#    r=[]
#    eff=[]
#    files = [f for f in os.listdir('.') if os.path.isfile(f)]
#    for file in files:
#        line = file.split("_")
#        if line[0] == 'Ngamma':
#            print line[1], getEff( line[1] )
#            r.append( line[1] )
#            eff.append( getEff( line[1] ) )
#            
#    plt.plot(r,eff,'*')
#    plt.show()
#
#plot_eff()
#
##id=1
##ein = read_data( 'Injection_'+str(id) )
##eout = read_data( 'Ngamma_'+str(id) )
##ein[1] *= ein[0]*ein[0]
##eout[1] *= eout[0]*eout[0]
##plt.loglog(ein[0],ein[1])
##plt.loglog(eout[0],eout[1])
##plt.show()
