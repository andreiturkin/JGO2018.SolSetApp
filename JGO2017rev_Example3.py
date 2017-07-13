from numpy import sin, cos
# Pi
from math import pi
#Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from JGO2017rev_NUC import CoveringTree, Box
from JGO2017rev_Plotting import PlottingTree

class Example3(CoveringTree, PlottingTree):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta=0, ShowCovPrc=False, bIntl=True):
        # Define the Initial Rectangle P
        side = 2*pi
        iBox = Box((pi/2.0, -pi), (side, side))
        pBox = Box((pi/2.0, -pi), (side, side))

        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, ShowCovPrc)

    def getResProcessedLevels(self):
        return super(Example3, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example3, self).getResIterations()

    def SaveResults(self, fileName, AddRings):
        super(Example3, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=AddRings)

    def SaveSolution(self, fName):
        super(Example3, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example3, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example3, self).LoadSolution(fName)

############################################################################################
# Private Methods
############################################################################################

    ########################################################################
    #                               Rho Constraints
    ########################################################################
    @staticmethod
    def phi(x):
        def g1(x):
            return sin(x[0])

        def g2(x):
            return -cos(x[1])

        return max(g1(x), g2(x))

    def getPhi(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        cphival = self.phi(((xmin+xmax)/2.0, (ymin+ymax)/2.0))
        L = 1.0
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
