#! /usr/bin/env python3

class blazarPlotter:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rc,rcParams
    import sys
    import math
    processList = "syn", "ssc", "ext1", "ext2", "accd", "blrPl", "blrSp", "blrGu", "dtPl", "dtSp", "dtGu", "accdGu", "QAccd"

    class processData:
        import numpy as np
        def __init__( self, _id ):
            self.name = _id
            self.xPeak = 0.0
            self.peakLuminosity = 0.0
            self.bolometricLuminosity = 0.0
            self.x = []
            self.y = []
            self.yNoDisk = []
            self.xVecPeak = []
            self.yVecPeak = []
            self.ffintLOG = []
            self.XLinIndex = 0.0
            self.XTwoPointIndex = 0.0

    def getPlt( self ):
        return self.plt

    def __init__( self, _directory, _fileType ):
        self.directory = _directory
        self.filetype = _fileType
        self.processFound = []
        self.yAxisLabel = r'$v L_{v}$'
        self.xAxisLabel = r'$v$'
        self.bolometricEtaRad = 0.0
        self.peakEtaRad = 0.0
        self.totalPeakLuminosity = 0.0
        self.totalBolometricLuminosity = 0.0

    def getBroadbandSpectra( self ):
        vmin = 1.0e100
        vmax = 1.0
        precision = 10000
        for process in self.processFound:
            if self.np.min( process.x ) < vmin:
                vmin = self.np.min( process.x )
            if self.np.max( process.x ) > vmax:
                vmax = self.np.max( process.x )
        x = self.np.logspace( self.np.log10(vmin), self.np.log10(vmax), precision )
        y = self.np.zeros( len( x ) )
        yNoDisk = self.np.zeros( len( x ) )

        from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
        ffit = []
        for process in self.processFound:
            process.ffintLOG = InterpolatedUnivariateSpline( self.np.log10(process.x), self.np.log10(process.y) )

        for i in range( len(x) ):
            for process in self.processFound:
                if x[i] >= self.np.min(process.x) and x[i] <= self.np.max(process.x):
                    if process.name != 'QAccd':
                        yNoDisk[i] += self.math.pow( 10.0, process.ffintLOG( self.np.log10(x[i]) ) )
                    y[i] += self.math.pow( 10.0, process.ffintLOG( self.np.log10(x[i]) ) )
        
        procObj = self.processData( 'obs' )
        self.processFound.append( procObj )
        x = self.movAvg( x, 5 )
        y = self.movAvg( y, 5 )
        yNoDisk = self.movAvg( yNoDisk, 5 )
        x = self.movAvg( x, 10 )
        y = self.movAvg( y, 10 )
        yNoDisk = self.movAvg( yNoDisk, 10 )
        x = self.movAvg( x, 20 )
        y = self.movAvg( y, 20 )
        yNoDisk = self.movAvg( yNoDisk, 20 )
        procObj.x = x
        procObj.y = y
        procObj.yNoDisk = yNoDisk
        self.plotData( procObj )
        self.getBolometricLuminosity( procObj )
        self.findPeaksBroadband( procObj )
        self.calcIndex( procObj, 4.84e17, 2.42e18 )

    def checkProcess( self ):
        import os.path
        for process in self.processList:
            if self.filetype == "LvPointAvg":
                self.filename = self.directory+"/"+self.filetype+"_"+process
            else:
                self.filename = self.directory+"/"+self.filetype+"_"+process+"_1"
            if os.path.isfile( self.filename ):
                procObj = self.processData( process )
                self.processFound.append( procObj )
                self.readData( procObj )
                self.plotData( procObj )
                self.findPeaks( procObj )
                self.plotPeaks( procObj )
                self.getBolometricLuminosity( procObj )

    def getEtaRadiation( self, Gamma, etaEle, etaDiss, Ljet0 ):
        for process in self.processFound:
            self.totalBolometricLuminosity += process.bolometricLuminosity
            self.totalPeakLuminosity += process.peakLuminosity

        self.peakEtaRad = self.totalPeakLuminosity/( self.math.pow(Gamma,2)*etaEle*etaDiss*Ljet0 )
        self.bolometricEtaRad = self.totalBolometricLuminosity/( self.math.pow(Gamma,2)*etaEle*etaDiss*Ljet0 )

    def getBolometricLuminosity( self, process ):
        integral=0
        for i in range(len(process.x)-1):
            integral += 0.5*(process.y[i]/process.x[i]+process.y[i+1]/process.x[i+1])*(process.x[i+1]-process.x[i])
        process.bolometricLuminosity = integral

    def readData( self, process ):
        file = open( self.filename, 'r' )
        _x = []
        _y = []
        for line in file.readlines():
            _x.append( float(line.split(" ")[0]) )
            if process.name == "QAccd":
                _y.append( float(line.split(" ")[1]) )
            else:
                _y.append( float(line.split(" ")[2]) )

        # throw away points with y=0
        for i in range( len(_y) ):
            if _y[i] > 0.0:
                process.y.append( _y[i] )
                process.x.append( _x[i] )

        process.x = self.np.array( process.x )
        process.y = self.np.array( process.x*process.y )
        process.x = self.movAvg( process.x, 2 )
        process.y = self.movAvg( process.y, 2 )
        process.x = self.movAvg( process.x, 4 )
        process.y = self.movAvg( process.y, 4 )
        file.close()

    def plotData( self, process ):
        self.plt.ylabel( self.yAxisLabel )
        self.plt.xlabel( self.xAxisLabel )
        self.plt.ylim( 1.0e-4*self.np.max( process.y ), 5.0*self.np.max( process.y ) )
        self.plt.loglog( process.x, process.y, label = process.name )

    def plotUpdate( self ):
        self.plt.legend( loc=2 )
        self.plt.show( )

    def movAvg( self, data, WINDOW ):
        import numpy as np
        weightings = np.repeat(1.0, WINDOW) / WINDOW
        movavg = np.convolve(data, weightings)[WINDOW-1:-(WINDOW-1)]
        return movavg

    def findPeaks( self, process ):
        for i in range( len(process.x)-2 ):
            if( process.y[i+1] > process.y[i] and process.y[i+2] < process.y[i+1] ):
                process.xPeak = process.x[i+1]
                process.peakLuminosity = process.y[i+1]
                break

    def findPeaksBroadband( self, process ):
        lowerPeakBorder = 5.0e14
        upperPeakBorder = 1.0e16
        for i in range( len(process.x)-2 ):
            if( process.yNoDisk[i+1] > process.yNoDisk[i] and process.yNoDisk[i+2] < process.yNoDisk[i+1] ):
                process.xVecPeak.append( process.x[i+1] )
                process.yVecPeak.append( process.yNoDisk[i+1] )
                if len( process.xVecPeak ) == 2:
                    break

    def plotPeaks( self, process ):
        self.plt.loglog( process.xPeak, process.peakLuminosity, 'o')

    def calcIndex( self, process, x1, x2 ):

        x_low = float( x1 )
        x_high = float( x2 )
        
        x_fit = []
        y_fit = []
        for i in range( len(process.x) ):
            if x_low <= process.x[i] and process.x[i] <= x_high:
                x_fit.append( self.np.log10( process.x[i] ) )
                y_fit.append( self.np.log10( process.y[i] ) )

        if len( x_fit ) < 10:
            pass
        else:
            [ [a, b], CM ] = self.np.polyfit( x_fit, y_fit, 1, full=False, cov=True )
            process.XLinIndex = 1.0-a
            process.XTwoPointIndex = 1.0-( (y_fit[-1]-y_fit[0])/(x_fit[-1]-x_fit[0]) )

#    db = self.math.sqrt( self.np.diag( CM )[1] )
#
#         def f( x, a, b ):
#             y = []
#             for i in x:
#                 y.append( a*i+b )
#             return y
#
#    y_fitted = f( x_fit, a, b )
#
#    #return [-(a-1.0),da,self.np.power(10,x_fit),self.np.power(10,y_fitted)]


# ------------------------------------------------------------------------
import sys

dirName = sys.argv[1]
fileType = sys.argv[2]
ifPlot = sys.argv[3]

graph = blazarPlotter( dirName, fileType )
graph.checkProcess( )
# Gamma, etaEle, etaDiss, Ljet0
graph.getEtaRadiation( 15.0, 0.5, 0.1, 3.25e+46 )
graph.getBroadbandSpectra( )

for process in graph.processFound:
    if process.name == 'obs':
        print( " %s %le %le %le" % (process.name, process.bolometricLuminosity, process.XLinIndex, process.XTwoPointIndex), end='' )
        for i in range( len(process.xVecPeak) ):
            print( " %le %le" % (process.xVecPeak[i], process.yVecPeak[i]), end='' )
    else:
        print( " %s %le %le %le" % (process.name, process.xPeak, process.peakLuminosity, process.bolometricLuminosity), end='' )


print( " %le %le %le %le" % (graph.peakEtaRad, graph.bolometricEtaRad, graph.totalPeakLuminosity, graph.totalBolometricLuminosity) )

if ifPlot == 'True':
    p = graph.getPlt( )
    graph.plotUpdate( )    

