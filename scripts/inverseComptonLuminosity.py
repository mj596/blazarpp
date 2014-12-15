#! /usr/bin/env python

class electrons:
    import numpy as np
    import math
    import matplotlib.pyplot as plt

    def __init__( self ):
        self.p1 = 0.5
        self.p2 = 2.5
        self.gammabreak = 100.0
        self.gammamin = 1.0
        self.gammamax = 1.0e4
        self.K = 1.0e49
    
    def info( self ):
        print( 'Initializing electrons.' )
        print( 'p1 = %1.1f' % self.p1 )
        print( 'p2 = %1.1f' % self.p2 )
        print( 'gammabreak = %.1le' % self.gammabreak )
        print( 'gammamin = %.1le' % self.gammamin )
        print( 'gammamax = %.1le' % self.gammamax )
        print( 'K = %.1le' % self.K )

    def getgamma( self ):
        return self.gamma

    def getNgamma( self ):
        return self.Ngamma

    def setK( self,_k ):
        self.K = _k

    def setp1( self,_p1 ):
        self.p1 = _p1

    def setp1( self, _p2 ):
        self.p2 = _p2

    def setgammabreak( self, _gb ):
        self.gammabreak = _gb

    def setgammamin( self, _gmin ):
        self.gammamin = _gmin

    def setgammamax( self, _gmax ):
        self.gammamax = _gmax

    def fQ( self, g ):
        if g >= self.gammamin and g <= self.gammamax:
            if g <= self.gammabreak:
                return self.K*self.math.pow( self.gammabreak, self.p1 - self.p2 )*self.math.pow( g, -self.p1 )
            if g > self.gammabreak:
                return self.K*self.math.pow( g, -self.p2 )

    def init( self ):
        self.gamma = self.np.logspace( self.math.log10(self.gammamin), self.math.log10(self.gammamax), 100 )
        self.Ngamma = []
        for i in self.gamma:
            self.Ngamma.append( self.fQ( i ) )
            
        self.plt.loglog( self.gamma, self.np.array(self.gamma)*self.np.array(self.gamma)*self.np.array(self.Ngamma) )
#        self.plt.show( )

class external:
    import numpy as np
    import math
    import matplotlib.pyplot as plt

    eV2erg = 1.60217646e-12
    mec2 = 8.18726233e-7
    PLANCK_H = 6.626e-27

    def __init__( self ):
        self.e = 10.0
        self.up = 1.0
    
    def info( self ):
        print( 'Initializing external radiation.' )
        print( 'e = %1.1f' % self.e )
        print( 'up = %0.1le' % self.up )

    def sete( self,_e ):
        self.e = _e

    def setup( self,_up ):
        self.up = _up

    def gete( self ):
        return self.e

    def getvext( self ):
        return self.e*self.eV2erg/self.PLANCK_H

    def getup( self ):
        return self.up

class calculator:
    import numpy as np
    import math
    import matplotlib.pyplot as plt

    eV2erg = 1.60217646e-12
    mec2 = 8.18726233e-7
    PLANCK_H = 6.626e-27
    SIGMA_T = 6.652453e-25

    theta = 0.0667
    Gamma = 15.0

    def __init__( self ):
        print( 'Initializing calculator.' )
        self.miu = -( self.math.cos(self.theta)-self.beta(self.Gamma) )/( 1.0-self.beta(self.Gamma)*self.math.cos(self.theta) )

    def ff( self, g, e0, e, miu ):
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

    def beta( self, _g ):
        return self.math.sqrt( self.math.pow(_g,2) - 1.0 )/_g

    def calculateICLuminosity( self, v, ele, ext ):
        e0 = ext.getvext( )*self.PLANCK_H*self.Gamma/self.mec2

        # integral
        Int = 0.0
        for i in range( len(ele.getgamma( )) ):
            e = self.PLANCK_H*v/self.mec2
            g = ele.getgamma( )[i]
            print( "%le %le" % (ele.getgamma()[i], ele.getNgamma( )[i]) )
            Int += ele.getNgamma( )[i]/ele.getgamma( )[i]*self.ff( g, e0, e, self.miu )
        print( 'Int %.1le' % Int )
        
#            Int += ele.getNgamma( )[i]/ele.getgamma( )[i]

        Int *= ele.getgamma( )[1]/ele.getgamma( )[0]
        Int *= ext.getup( )/self.math.pow( ext.getvext(), 2 )
        Int *= 3.0*v*self.SIGMA_T/(4.0*self.math.pow( self.Gamma, 2 ) )

# --------------------------------------------------------------------
    
ele = electrons( )
ele.info( )
ele.init( )
ext = external( )
ext.info( )
calc = calculator( )
calc.calculateICLuminosity( 1.0489e+20, ele, ext )

##define mec2  (8.18726233e-7)
#
#
#eV2erg = 1.60217646e-12
#mec2 = 8.18726233e-7
#PLANCK_H = 6.626e-27
#vext = 10.0*eV2erg/PLANCK_H
#
#gamma = np.logspace(2,2,1)
#v = np.logspace(15,30,1000)
#theta = 0.0667
#Gamma = 15.0
#
#e0 = vext*PLANCK_H*Gamma/mec2
#miu = -( math.cos(theta)-beta(Gamma) )/( 1.0-beta(Gamma)*math.cos(theta) )
#
#
#import matplotlib.pyplot as plt
#for g in gamma:
#    print g
#    fvec = []
#    for i in v:
#        e = PLANCK_H*i/mec2
#        fvec.append( f(g,e0,e,miu) )
#    plt.loglog( v, fvec )
#plt.show()
