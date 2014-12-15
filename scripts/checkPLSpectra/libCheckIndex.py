#!/usr/bin/env python

class checkIndex:
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import math

    def __init__(self, _plotType, _filename, _canvas):
        if _plotType == 'ele' or _plotType == 'lum':
            self.plotType = _plotType
            self.filename = _filename
            self.readData()
            self.canvas = _canvas
            self.setRcParams()
            self.BINSIZE = 20
            self.R2THRESHOLD = 0.95

        else:
            print "[libCheckIndex] Please use \'ele\' or \'lum\' instead of", _plotType,"-> Quit."
            self.sys.exit()

    def setRcParams(self):
	self.plt.rcParams.update({'legend.labelspacing':0.05})
	self.plt.rcParams.update({'legend.fontsize': 14})
	self.plt.rcParams.update({'legend.fontfamily':'monospace'})
	self.plt.rcParams.update({'font.size':14})
	self.plt.rcParams.update({'font.family':'monospace'})
	self.plt.rcParams.update({'axes.fontfamily':'monospace'})
	self.plt.rcParams.update({'axes.fontsize':14})
    
    def readData(self):
        self.x=[]
        self.y=[]
        file = open(self.filename,'r')

        _x=[]
        _y=[]

        if self.plotType == 'ele':
            for line in file.readlines():
                _x.append(float(line.split(" ")[0]))
                _y.append(float(line.split(" ")[1]))

            # throw away points with y=0
            for i in range(len(_y)):
                if _y[i] > 0.0:
                    self.y.append(_y[i])
                    self.x.append(_x[i])

            self.x = self.np.array(self.x)
            self.y = self.np.array(self.y)

            self.y *= self.x*self.x
            self.plotLabel = 'Blazar++ ele spectrum'
            self.xAxisLabel = '\$\gamma^2 N_{\gamma}\$'
            self.yAxisLabel = '\$\gamma}\$'

        if self.plotType == "lum":
            for line in file.readlines():
                _x.append(float(line.split(" ")[0]))
                _y.append(float(line.split(" ")[2]))

            # throw away points with y=0
            for i in range(len(_y)):
                if _y[i] > 0.0:
                    self.y.append(_y[i])
                    self.x.append(_x[i])

            self.x = self.np.array(self.x)
            self.y = self.np.array(self.y)

            self.y *= self.x
            self.plotLabel = 'Blazar++ lum spectrum'
            self.xAxisLabel = 'frequency'
            self.yAxisLabel = 'lum'
            
        file.close()

    def plot(self):
    	self.plt.ylabel(self.yAxisLabel)
    	self.plt.xlabel(self.xAxisLabel)
    	self.plt.title(self.plotLabel)
#        if self.plotType == 'ele':
#            self.plt.xlim(0.5*self.x[0],1.1*self.np.max(self.x))
#            self.plt.ylim(0.5*self.y[0],1.5*self.np.max(self.y))
#        if self.plotType == 'lum':
#            self.plt.ylim(0.5*self.y[0],5.0*self.np.max(self.y))
        self.plt.loglog(self.x,self.y,linewidth=2)
        
    def s2alpha(self,s):
        return (float(s)-1.0)/2.0

    def alpha2s(self,alpha):
        return 2.0*float(alpha)+1.0

    def addIndex(self,_array,_indexType):
        if len(_array) == 0:
            print "[libCheckIndex] No index specified! Quit."
            self.sys.exit()
        else:
            if _indexType == 'ele' or _indexType == 'lum':
                self.INDEX = []
                self.indexType = _indexType

                print "Adding index:",
                for i in _array:
                    if self.plotType == 'ele' and self.indexType == 'ele':
                        self.INDEX.append( i )
                        print i,
                    if self.plotType == 'lum' and self.indexType == 'lum':
                        self.INDEX.append( i )
                        print i,
                    if self.plotType == 'ele' and self.indexType == 'lum':
                        self.INDEX.append( self.alpha2s(i) )
                        print self.alpha2s(i),'(',i,')',
                    if self.plotType == 'lum' and self.indexType == 'ele':
                        self.INDEX.append( self.s2alpha(i) )
                        print self.s2alpha(i),'(',i,')',
                    
                self.INDEX = self.np.array(self.INDEX)
                print ""

            else:
                print "[libCheckIndex] Please use \'ele\' or \'lum\' instead of", _indexType, "-> Quit."
                self.sys.exit()

    def lineplot(self,x,_index,xp):
        def find_point(x,xp):
            from scipy import interpolate
            f = interpolate.interp1d(x[0], x[1])
            return f(xp)

        vy=[]
        
        if self.plotType == 'ele':
            index = float(_index)-2.0

        if self.plotType == "lum":
            index = float(_index)-3.0

        for i in x[0]:
            vy.append(10**(-index*(self.math.log10(i)-self.math.log10(xp))+self.math.log10(find_point(x,xp))))
        return self.np.array(vy)

    def getR2(self,x,index):
        # get mean value
        xmean = self.np.mean(x[0])
        ymean = self.np.mean(x[1])
        # get model values
        ymodel = self.lineplot(x,index,xmean)
        # caluclate R2
        y=x[1]
        r2up=0.0
        r2down=0.0
        for i in range(len(y)):
            r2up += self.math.pow(y[i]-ymodel[i],2.0)
            r2down += self.math.pow(y[i]-ymean,2.0)
    
        return 1.0-(r2up/r2down)

    def findIndex(self):
        for j in self.INDEX:
            print "Checking index", j, '->',
            INDEX_FLAG = True
            xBIN = []
            yBIN = []
            for i in range( len(self.x)/self.BINSIZE ):
                if INDEX_FLAG:
                    xBIN.append( self.x[i*self.BINSIZE:(i+1)*self.BINSIZE-1] )
                    yBIN.append( self.y[i*self.BINSIZE:(i+1)*self.BINSIZE-1] )
#                    print xBIN[i], yBIN[i], self.getR2([xBIN[i],yBIN[i]],j), j # DEBUG
#                    self.plt.loglog(xBIN[i],yBIN[i],linewidth='6') # DEBUG
#                    self.plt.loglog(self.x,self.lineplot([self.x,self.y],j,self.np.mean(xBIN[i]))) # DEBUG
                    if self.getR2([xBIN[i],yBIN[i]],j) > self.R2THRESHOLD:
                        print 'Found:', 'range',(xBIN[i])[0],'-',(xBIN[i])[-1],'R2:',self.getR2([xBIN[i],yBIN[i]],j)
                        self.plt.loglog(xBIN[i],yBIN[i],linewidth='6')
                        if self.plotType == 'ele':
                            _label = 's='+str(j)+' (alpha='+str(self.s2alpha(j))+')'
                        if self.plotType == 'lum':
                            _label = 'alpha='+str(j)+' (s='+str(self.alpha2s(j))+')'
                        self.plt.loglog(self.x,self.lineplot([self.x,self.y],j,self.np.mean(xBIN[i])),label=_label)
                        INDEX_FLAG = False
            if INDEX_FLAG:
                print "Not found."
        self.plt.legend(loc=8)
