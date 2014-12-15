#! /usr/bin/env python3
class plotUpeR( ):
    import os
    import sys
    import matplotlib.pyplot as plt
    from matplotlib import rc,rcParams
    processList = "syn", "ssc", "blrSp", "dtSp", "ext1", "ext2", "accd", "blrPl", "dtPl", "blrGu", "dtGu", "accdGu"

    def __init__( self, _directory ):
        self.directory = _directory + "/out"
        print( "-> Using %s" % (self.directory) )
        self.processFound = []

    def checkRunList( self ):
        self.runList = self.os.listdir( self.directory )
        print( "-> Found %d runs" % len( self.runList ) )
        
    def checkProcess( self ):
        for process in self.processList:
            self.filename = self.directory + "/" + self.runList[0] + "/UpeR_" + process
            if self.os.path.isfile( self.filename ):
                print ( "-> Found %s." % process )
                self.processFound.append( process )

    def readData( self ):
        self.x = []
        self.y = []
        file = open( self.filename, 'r' )
        _x=[]
        _y=[]
        for line in file.readlines():
            _x.append( float(line.split(" ")[0]) )
            _y.append( float(line.split(" ")[1]) )

        # throw away points with y=0
        for i in range( len(_y) ):
            if _y[i] > 0.0:
                self.y.append( _y[i] )
                self.x.append( _x[i] )

        file.close()

    def plot( self ):
        for process in self.processFound:
            print( "-> Plotting %s" % (process) )
            self.initializeDataVector( )
            for run in self.runList:
                self.filename = self.directory + "/" + run + "/UpeR_" + process
                self.readData( )
                self.updateDataVector( )
            self.sortSumData( )
            self.plotSumData( process )
        self.plotUpdate( )

    def plotData( self, _process ):
        self.plt.loglog( self.x, self.y )

    def plotSumData( self, _process ):
        self.plt.loglog( self.x_sum, self.y_sum, '.', label = _process )

    def plotUpdate( self ):
        self.plt.legend( loc=2 )
        self.plt.show( )

    def initializeDataVector( self ):
        self.x_sum = []
        self.y_sum = []

    def updateDataVector( self ):
        for i in range( len(self.x) ):
            self.x_sum.append( self.x[i] )
            self.y_sum.append( self.y[i] )

    def sortSumData( self ):
        print( "No sorting yet!" )
        return 0


import sys
a = plotUpeR( sys.argv[1] )
a.checkRunList( )
a.checkProcess( )
a.plot( )
    
    
