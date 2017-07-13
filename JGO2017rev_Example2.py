# Math functions and constants
from math import sqrt, sin, cos, acos, pi
# Interval analysis
from interval import interval
# Arrays
import numpy as np

from JGO2017rev_NUC import CoveringTree, Box
from JGO2017rev_Plotting import PlottingTree
# from JGO2017rev_Example3_Phi import Phi

########################################################################
#   Table of robot parameters *
#   i               1               2               3
#-----------------------------------------------------------------------
#   x_ai            -15             15              0
#   y_ai            -5sqrt(3)       -5sqrt(3)       10sqrt(3)
#   x_bi            -5              5               0
#   y_bi            -5sqrt(3)/3     -5sqrt(3)/3     10sqrt(3)/3
#   rho_min         12              12              12
#   rho_max         27              27              27
#-----------------------------------------------------------------------
#   * see Gosselin C., Jean M. - Determination of the workspace of planar
#   parallel manipulators with joint limits.
########################################################################
# For the constant angle:
# angle = 50
#

# class Example2(CoveringTree, PlottingTree, Phi):
class Example2(CoveringTree, PlottingTree):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta=0, iangle=(10.0/180.0)*pi, ShowCovPrc=False):
        # Full Range
        # Define the Initial Rectangle P
        corner = -20.0
        side = 40.0
        iBox = Box((corner, corner), (side, side))
        pBox = Box((corner, corner), (side, side))

        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, iangle, ShowCovPrc)
        # Phi.__init__(self)

        # Initialize the parameters of the planar parallel robot
        # A constraints
        self.__a = []
        self.__a.append(np.array([-15.0, -5.0*sqrt(3.0)]))
        self.__a.append(np.array([15.0, -5.0*sqrt(3.0)]))
        self.__a.append(np.array([0.0, 10.0*sqrt(3.0)]))
        # B constraints
        self.__b = []
        self.__b.append(np.array([-5.0, -5.0*sqrt(3.0)/3.0]))
        self.__b.append(np.array([5.0, -5.0*sqrt(3.0)/3.0]))
        self.__b.append(np.array([0.0, 10.0*sqrt(3.0)/3.0]))

        self.__rho = [12.0, 27.0]
        self.__angle = iangle
        self.__p = [self.__a[u][0] - self.__b[u][0]*cos(self.__angle) + self.__b[u][1]*sin(self.__angle) for u in [0,1,2]]
        self.__q = [self.__a[u][1] - self.__b[u][0]*sin(self.__angle) - self.__b[u][1]*cos(self.__angle) for u in [0,1,2]]

    def getResProcessedLevels(self):
        return super(Example2, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example2, self).getResIterations()

    def SaveResults(self, fileName, AddRings):
        super(Example2, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=AddRings)

    def SaveSolution(self, fName):
        super(Example2, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example2, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example2, self).LoadSolution(fName)

############################################################################################
# Private Methods
############################################################################################
    def getRho(self):
        return self.__rho

    def getCi(self, angle, i):
        return self.__a[i] - self.RotL(self.__b[i], angle)

    # Rotate u counterclockwise
    @staticmethod
    def RotL(u, a):
        return np.array([u[0]*cos(a) - u[1]*sin(a), u[0]*sin(a) + u[1]*cos(a)])

    ########################################################################
    #                               Rho Constraints
    ########################################################################
    def gRhoi(self, x, angle, i):
        u_c1 = self.__a[i]- self.RotL(self.__b[i], angle)
        rho2i = np.linalg.norm(x-u_c1)**2
        return np.array([rho2i - self.__rho[1]**2, self.__rho[0]**2 - rho2i])

############################################################################################
# Abstract Methods
############################################################################################

    ########################################################################################
    # CoveringTree
    ########################################################################################
    def getMinMaxVal(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        x = [(xmin+xmax)/2.0, (ymin+ymax)/2.0]
        flist = ['gRhoi']

        res = [r for f in flist for i in [0, 1, 2] for r in getattr(self, f)(x, self.__angle, i)]
        phi = max(res)
        # Lc = self.CasADi_getLipConst(iBox)
        Lg = [2.0*sqrt(max((xmin-self.__p[u])**2, (xmax-self.__p[u])**2) + \
                       max((ymin-self.__q[u])**2, (ymax-self.__q[u])**2)) \
                       for u in [0, 1, 2]]
        L = max(Lg)
        # L = Lc
        phiMin = phi - (float(L)/2.0)*iBox.getDiam()
        phiMax = phi + (float(L)/2.0)*iBox.getDiam()
        # return minimum and maximum values
        return phiMin, phiMax

    ########################################################################################
    # PlottingTree
    ########################################################################################
    def AdditionalPlotting(self, ax):
        ci = [self.getCi(self.vals[0], i) for i in [0, 1, 2]]
        for c in ci:
            for r in self.getRho():
                self.drawCircle(c, r)
        return

############################################################################################
# !Abstract Methods
############################################################################################
