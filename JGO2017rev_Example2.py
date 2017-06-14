import numpy as np
from numpy import linalg as la
import datetime
# Pi
from math import pi
# Interval analysis
from interval import interval
#Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#Saving data to files
import os
import json

from JGO2017rev_NUC import CoveringTree, Box
from JGO2017rev_Plotting import PlottingTree
from JGO2017rev_Example2_Phi import Phi

class Example2(CoveringTree, PlottingTree, Phi):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta=0, ShowCovPrc=False):
        # Zoom up the bottom triangle
        # iBox = Box((-8, -18, (8.0/180.0)*pi), (14, 8, (12.0/180.0)*pi))

        # Full Range
        # Define the Initial Rectangle P
        corner = -20.0
        side = 40.0
        iBox = Box((corner, corner, (9.0/180.0)*pi), (side, side, (2.0/180.0)*pi))
        pBox = Box((corner, corner, (9.0/180.0)*pi), (side, side, (2.0/180.0)*pi))

        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, ShowCovPrc)
        Phi.__init__(self)

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
        ci = [self.getCi(self.vals[0], i) for i in [0, 1, 2]]

        for c in ci:
            for r in self.getRho():
                self.drawCircle(c, r)

        theta_a_bounds = self.getTheta_a_Bounds()
        theta_b_bounds = self.getTheta_b_Bounds()
        theta_a_intl = [interval[theta_a_bounds[0][i], theta_a_bounds[1][i]] for i in [0, 1, 2]]
        theta_b_intl = [interval[theta_b_bounds[0][i] + self.vals[0] - pi,\
                                 theta_b_bounds[1][i] + self.vals[0] - pi] for i in [0, 1, 2]]

        psi = [(theta_b_intl[0][0][0], theta_b_intl[0][0][1]),\
               (theta_b_intl[1][0][0], theta_b_intl[1][0][1]),\
               (theta_b_intl[2][0][0], theta_b_intl[2][0][1])]

        for idx, c in enumerate(ci):
            for r in self.getRho():
                self.drawArc(c, r, psi[idx])

        for idx, c in enumerate(ci):
            self.drawLine(c, self.getRho(), psi[idx])
        return

############################################################################################
# !Abstract Methods
############################################################################################
