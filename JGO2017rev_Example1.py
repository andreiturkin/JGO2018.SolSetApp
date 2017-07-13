from math import sqrt

from JGO2017rev_NUC import Box, CoveringTree
from JGO2017rev_Plotting import PlottingTree
from JGO2017rev_GlobOpt import GlobalOptimization

def phi(x):
    def g2(x):
        return x[0]**2.0 + x[1]**2.0 - 1

    def g1(x):
        return 0.999**2 - x[0]**2.0 - x[1]**2.0
    return max(g1(x), g2(x))

def getLipVal(iBox):
    bounds = iBox.getBounds()
    xmin = bounds[0][0]
    xmax = bounds[0][1]
    ymin = bounds[1][0]
    ymax = bounds[1][1]

    return 2*sqrt(max(abs(xmin), abs(xmax))**2 + max(abs(ymin), abs(ymax))**2)

class Example1_GlobOpt(CoveringTree, PlottingTree):
    def __init__(self, idelta=0, ShowCovPrc=False):
        # The Initial Rectangle P definition
        corner = -1.5
        side = 3.0

        # The algorithm parameters initialization
        self.delta = idelta
        self.eps = (getLipVal(Box((corner, corner), (side, side)))*idelta)/2.0

        # The CoveringTree class constructor
        CoveringTree.__init__(self, Box((corner, corner), (side, side)), self.delta, ShowCovPrc, self.eps)
        PlottingTree.__init__(self, Box((corner, corner), (side, side)), ShowCovPrc=ShowCovPrc)

    def SaveResults(self, fileName, AddRings):
        super(Example1_GlobOpt, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=AddRings)

    def SaveSolution(self, fName):
        super(Example1_GlobOpt, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example1_GlobOpt, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example1_GlobOpt, self).LoadSolution(fName)

    def getResProcessedLevels(self):
        return super(Example1_GlobOpt, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example1_GlobOpt, self).getResIterations()

    def AdditionalPlotting(self, ax):
        self.drawCircle((0,0), 1)
        self.drawCircle((0,0), 0.999)

    def getMinMaxVal(self, iBox):
        gopt = GlobalOptimization(phi, getLipVal, self.eps)
        return gopt.getMinVal(iBox),\
               gopt.getMaxVal(iBox)

class Example1_AppxGlobL(CoveringTree, PlottingTree):
    def __init__(self, idelta=0, ShowCovPrc=False):
        #Define the Initial Rectangle P
        corner = -1.5
        side = 3.0
        self.__L = getLipVal(Box((corner, corner), (side, side)))
        CoveringTree.__init__(self, Box((corner, corner), (side, side)), idelta)
        PlottingTree.__init__(self, Box((corner, corner), (side, side)), ShowCovPrc=ShowCovPrc)

    def getResProcessedLevels(self):
        return super(Example1_AppxGlobL, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example1_AppxGlobL, self).getResIterations()

    def SaveResults(self, fileName, AddRings):
        super(Example1_AppxGlobL, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=AddRings)

    def SaveSolution(self, fName):
        super(Example1_AppxGlobL, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example1_AppxGlobL, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example1_AppxGlobL, self).LoadSolution(fName)

    def AdditionalPlotting(self, ax):
        #Example 1:
        self.drawCircle((0,0), 1)
        self.drawCircle((0,0), 0.999)

    def getMinMaxVal(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        phival = phi(((xmin+xmax)/2.0, (ymin+ymax)/2.0))
        return phival - (self.__L/2.0)*iBox.getDiam(), \
               phival + (self.__L/2.0)*iBox.getDiam()

class Example1_AppxLocL(CoveringTree, PlottingTree):

    def __init__(self, idelta=0, ShowCovPrc=False):
        #Define the Initial Rectangle P
        corner = -1.5
        side = 3.0

        CoveringTree.__init__(self, Box((corner, corner), (side, side)), idelta)
        PlottingTree.__init__(self, Box((corner, corner), (side, side)), ShowCovPrc=ShowCovPrc)

    def getResProcessedLevels(self):
        return super(Example1_AppxLocL, self).getResProcessedLevels()

    def getResIterations(self):
        return super(Example1_AppxLocL, self).getResIterations()

    def SaveResults(self, fileName, AddRings):
        super(Example1_AppxLocL, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=AddRings)

    def SaveSolution(self, fName):
        super(Example1_AppxLocL, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(Example1_AppxLocL, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(Example1_AppxLocL, self).LoadSolution(fName)

    def AdditionalPlotting(self, ax):
        self.drawCircle((0,0), 1)
        self.drawCircle((0,0), 0.999)

    def getMinMaxVal(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        L = getLipVal(iBox)
        phival = phi(((xmin+xmax)/2.0, (ymin+ymax)/2.0))
        return phival - (L/2.0)*iBox.getDiam(), \
               phival + (L/2.0)*iBox.getDiam()
