from numpy import sin, cos, exp
# Pi
from math import pi
#Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from JGO2017rev_NUC import CoveringTree, Box
from JGO2017rev_Plotting import PlottingTree

class Example2_Lipz(CoveringTree, PlottingTree):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta=0, ShowCovPrc=False, bIntl=True):
        # Define the Initial Rectangle P
        side = pi
        iBox = Box((pi, 1.0), (side, side))
        pBox = Box((pi, 1.0), (side, side))

        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, ShowCovPrc)

    def getResProcessedLevels(self):
        return super(Example2_Lipz, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example2_Lipz, self).getResIterations()

    def SaveResults(self, fileName, Zoom=False, ZoomBox=None):
        super(Example2_Lipz, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=False, Zoom=Zoom, ZoomBox=ZoomBox)

    def SaveSolution(self, fName):
        super(Example2_Lipz, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example2_Lipz, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example2_Lipz, self).LoadSolution(fName)

############################################################################################
# Private Methods
############################################################################################

    ########################################################################
    #                               Rho Constraints
    ########################################################################
    @staticmethod
    def phi(x):
        def g1(x):
            return x[0]*sin(x[0]) + 0.1*(x[0]**2) + 1.0

        def g2(x):
            return cos(x[1]) + 0.1*(x[1]**2)

        return max(g1(x), g2(x))
    @staticmethod
    def maxnormnablaphi(x):
        def g1(x):
            return abs(imath.sin(x[0]) + x[0]*imath.cos(x[0]) + 0.2*x[0])

        def g2(x):
            return abs(-1.0*imath.sin(x[1]) + 0.2*x[1])

        res = [g1(x), g2(x)]
        normnablaphi = interval[max([interval.hull([interval(u) for u in x])[0][0] for x in res]),\
                                max([interval.hull([interval(u) for u in x])[0][1] for x in res])]
        return normnablaphi[0][1]

    def getPhi(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        cphival = self.phi(((xmin+xmax)/2.0, (ymin+ymax)/2.0))
        iBInt = [interval([xmin,xmax]),interval([ymin,ymax])]
        L = self.maxnormnablaphi(iBInt)
        phiMin = cphival - (L/2.0)*iBox.getDiam()
        phiMax = cphival + (L/2.0)*iBox.getDiam()

        return phiMin, phiMax

############################################################################################
# Abstract Methods
############################################################################################

    ########################################################################################
    # CoveringTree
    ########################################################################################
    def getMinMaxVal(self, iBox):

        phiMin, phiMax = self.getPhi(iBox)
        # return minimum and maximum values
        return phiMin, phiMax

    ########################################################################################
    # PlottingTree
    ########################################################################################
    def AdditionalPlotting(self, ax):
        return

############################################################################################
# !Abstract Methods
############################################################################################
import numpy as np
import itertools
import random
class Example2_LipzHcs(CoveringTree, PlottingTree):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta=0, ShowCovPrc=False, bIntl=True):
        # Define the Initial Rectangle P
        side = pi
        iBox = Box((pi, 1.0), (side, side))
        pBox = Box((pi, 1.0), (side, side))
        self.__initBoxDiam = iBox.getDiam()
        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, ShowCovPrc)

    def getResProcessedLevels(self):
        return super(Example2_LipzHcs, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example2_LipzHcs, self).getResIterations()

    def SaveResults(self, fileName, Zoom=False, ZoomBox=None):
        super(Example2_LipzHcs, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=False, Zoom=Zoom, ZoomBox=ZoomBox)

    def SaveSolution(self, fName):
        super(Example2_LipzHcs, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example2_LipzHcs, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example2_LipzHcs, self).LoadSolution(fName)

############################################################################################
# Private Methods
############################################################################################

    ########################################################################
    #                               Rho Constraints
    ########################################################################
    @staticmethod
    def phi(x):
        def g1(x):
            return x[0]*sin(x[0]) + 0.1*(x[0]**2) + 1.0

        def g2(x):
            return cos(x[1]) + 0.1*(x[1]**2)
        return max(g1(x), g2(x))

    def maxnormnablaphi(self, xmin, xmax, ymin, ymax):
        N = 10
        xpoints = [random.uniform(xmin, xmax) for _ in range(N)]
        ypoints = [random.uniform(ymin, ymax) for _ in range(N)]
        points = zip(xpoints, ypoints)

        maxval = 0

        retval = max([abs(self.phi(point1) - self.phi(point2))/np.linalg.norm(np.array(point1)-np.array(point2)) \
                      for idx1, point1 in enumerate(points) for idx2, point2 in enumerate(points) \
                      if idx1 != idx2])
        return retval

    def getPhi(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        cphival = self.phi(((xmin+xmax)/2.0, (ymin+ymax)/2.0))
        alpha = 1.2
        L = alpha*self.maxnormnablaphi(xmin, xmax, ymin, ymax)
        phiMin = cphival - (L/2.0)*iBox.getDiam()
        phiMax = cphival + (L/2.0)*iBox.getDiam()

        return phiMin, phiMax

############################################################################################
# Abstract Methods
############################################################################################

    ########################################################################################
    # CoveringTree
    ########################################################################################
    def getMinMaxVal(self, iBox):

        phiMin, phiMax = self.getPhi(iBox)
        # return minimum and maximum values
        return phiMin, phiMax

    ########################################################################################
    # PlottingTree
    ########################################################################################
    def AdditionalPlotting(self, ax):
        return

############################################################################################
# !Abstract Methods
############################################################################################

from interval import interval, imath
class Example2_Intl(CoveringTree, PlottingTree):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta=0, ShowCovPrc=False, bIntl=True):
        # Define the Initial Rectangle P
        side = pi
        iBox = Box((pi, 1.0), (side, side))
        pBox = Box((pi, 1.0), (side, side))

        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, ShowCovPrc)

    def getResProcessedLevels(self):
        return super(Example2_Intl, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example2_Intl, self).getResIterations()

    def SaveResults(self, fileName, Zoom=False, ZoomBox=None):
        super(Example2_Intl, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=False, Zoom=Zoom, ZoomBox=ZoomBox)

    def SaveSolution(self, fName):
        super(Example2_Intl, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example2_Intl, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example2_Intl, self).LoadSolution(fName)

############################################################################################
# Private Methods
############################################################################################

    ########################################################################
    #                               Rho Constraints
    ########################################################################
    @staticmethod
    def phi(x):
        def g1(x):
            return x[0]*imath.sin(x[0]) + 0.1*(x[0]**2) + 1.0

        def g2(x):
            return imath.cos(x[1]) + 0.1*(x[1]**2)

        res = [g1(x), g2(x)]
        phi = interval[max(interval.hull([interval(u) for u in x])[0][0] for x in res),\
                       max(interval.hull([interval(u) for u in x])[0][1] for x in res)]

        return phi[0][0], phi[0][1]

    def getPhi(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        iBInt = [interval([xmin,xmax]),interval([ymin,ymax])]
        oPhiInt_min, oPhiInt_max  = self.phi(iBInt)
        return oPhiInt_min, oPhiInt_max

############################################################################################
# Abstract Methods
############################################################################################

    ########################################################################################
    # CoveringTree
    ########################################################################################
    def getMinMaxVal(self, iBox):

        phiMin, phiMax = self.getPhi(iBox)
        # return minimum and maximum values
        return phiMin, phiMax

    ########################################################################################
    # PlottingTree
    ########################################################################################
    def AdditionalPlotting(self, ax):
        return

############################################################################################
# !Abstract Methods
############################################################################################



