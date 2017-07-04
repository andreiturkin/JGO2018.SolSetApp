# Math functions and constants
from math import sqrt, sin, cos, acos, pi
# Numpy
import numpy as np
# Operator
import operator
# Min Float
import sys
# GlobalOptimization by using Brute Force
from scipy import optimize
# Interval analysis
from interval import interval
# Automatic Differentiation
import casadi

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
# angle = 10
#

class Phi(object):
    def __init__(self):
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
        self.__angle = (10.0/180.0)*pi

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

    def getIntl(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]
        return (interval[xmin, xmax], interval[ymin, ymax])

    def getPhi(self, iBox):
        intl = self.getIntl(iBox)
        inx = [interval(i).midpoint for i in intl]
        x = [i[0][0] for i in inx[:2]]
        flist = ['gRhoi']

        res = [r for f in flist for i in [0, 1, 2] for r in getattr(self, f)(x, self.__angle, i)]
        phi = max(res)
        L = self.CasADi_getLipConst(iBox)
        phiMin = phi - (float(L)/2.0)*iBox.getDiam()
        phiMax = phi + (float(L)/2.0)*iBox.getDiam()

        return phiMin, phiMax

    def CasADi_getLipConst(self, iBox):
        def Rotate(u, a):
            return np.array([u[0]*casadi.cos(a) - u[1]*casadi.sin(a), u[0]*casadi.sin(a) + u[1]*casadi.cos(a)])
        def NormSquared(u):
            return u[0]**2+u[1]**2
        def Norm(u):
            return casadi.sqrt(u[0]**2+u[1]**2)
        def casadi_Rhoi_init(xc, i): # x = (x,y,angle)
            ret = []

            u_c1 = self.__a[i]- Rotate(self.__b[i], self.__angle)
            rho2i = NormSquared(xc[:2]-u_c1)
            ret.append(rho2i - self.__rho[1]**2)
            ret.append(self.__rho[0]**2 - rho2i)

            return ret

        xc = casadi.SX.sym('xc',2)
        flist = ['casadi_Rhoi_init']
        SXfs = [item for sublist in [locals()[f](xc,i) for f in flist for i in [0, 1, 2]] for item in sublist]

        maxSXfs = casadi.fmax(SXfs[5], casadi.fmax(SXfs[4],\
        casadi.fmax(SXfs[3], casadi.fmax(SXfs[2], casadi.fmax(SXfs[1], SXfs[0])))))

        # Test for Taking Derivatives
        # abc = casadi.SX.sym('abc',3)
        # abcfunc = casadi.Function('abcfunc', [abc[0], abc[1], abc[2]],[abc[0]**2 + abc[1]**2 + abc[2]**2],['a','b','c'],['abcnorm2'])
        # dabcda = abcfunc.jacobian('a','abcnorm2')
        # dabcdb = abcfunc.jacobian('b','abcnorm2')
        # dabcdc = abcfunc.jacobian('c','abcnorm2')
        # norm_dabc = casadi.sqrt(dabcda(*dabcda.sx_in())[0]**2+dabcdb(*dabcdb.sx_in())[0]**2+dabcdc(*dabcdc.sx_in())[0]**2)

        phifunc = casadi.Function('phifunc', [xc[0], xc[1]], [maxSXfs], ['x','y'], ['phi'])
        phi_dx = phifunc.jacobian('x','phi')
        phi_dy = phifunc.jacobian('y','phi')
        dphi_norm = casadi.sqrt(phi_dx(*phi_dx.sx_in())[0]**2 + \
                                phi_dy(*phi_dy.sx_in())[0]**2)
        dphi_norm_func = casadi.Function('dphinormfunc', [xc[0], xc[1]], [dphi_norm], ['x','y'], ['normnablaphi'])
        dphi_norm_func_inv = casadi.Function('dphinorminvfunc', [xc[0], xc[1]], [-dphi_norm], ['x','y'], ['normnablaphiinv'])
        dphi_norm_func_inv_jac_x = dphi_norm_func_inv.jacobian('x', 'normnablaphiinv')
        dphi_norm_func_inv_jac_y = dphi_norm_func_inv.jacobian('y', 'normnablaphiinv')

        # Finding a max value for every norm in the list
        def iFunc(inx, *func):
            funcret = func[0].call([inx[0],inx[1]])
            return funcret[0]

        def idFuncInv(inx, *dfunc):
            ret_x = dfunc[1].call([inx[0],inx[1]])[0]
            ret_y = dfunc[2].call([inx[0],inx[1]])[0]
            return np.array([ret_x, ret_y])

        ret = []
        intl = self.getIntl(iBox)
        xmin, xmax = intl[0][0][0], intl[0][0][1]
        ymin, ymax = intl[1][0][0], intl[1][0][1]
        steplength = 10.0
        xstep = (xmax-xmin)/steplength
        ystep = (ymax-ymin)/steplength
        rranges = (slice(xmin, xmax, xstep), slice(ymin, ymax, ystep))
        maxbrute = optimize.brute(iFunc, rranges, args=[dphi_norm_func_inv], full_output=True, finish=None)
        bnds = ((xmin, xmax), (ymin, ymax))
        ret = optimize.minimize(iFunc, maxbrute[0], args=tuple([dphi_norm_func_inv, dphi_norm_func_inv_jac_x, dphi_norm_func_inv_jac_y]), \
                                jac=idFuncInv, method='L-BFGS-B', bounds=bnds)
        # Finding a max value for these value,
        # which is a Lipschitz constant for phi function
        return iFunc(ret.x, dphi_norm_func)

    def GlobOptExt(self, iBox):
        def getRes(inx):
            x = inx[0:2]
            angle = inx[2]
            flist = ['gRhoi']
            return [r for f in flist for i in [0, 1, 2] for r in getattr(self, f)(x, angle, i)]

        def PhiFunc(inx):
            res = getRes(inx)
            return max(res)

        def PhiFuncInv(inx):
            res = getRes(inx)
            return -max(res)

        intl = self.getIntl(iBox)
        xmin, xmax = intl[0][0][0], intl[0][0][1]
        ymin, ymax = intl[1][0][0], intl[1][0][1]
        steplength = 10.0
        xstep = (xmax-xmin)/steplength
        ystep = (ymax-ymin)/steplength
        rranges = (slice(xmin, xmax, xstep), slice(ymin, ymax, ystep))
        minbrute = optimize.brute(PhiFunc, rranges, full_output=True, finish=None)
        maxbrute = optimize.brute(PhiFuncInv, rranges, full_output=True, finish=None)

        return minbrute[1], PhiFunc(maxbrute[0]), getRes(minbrute[0]), getRes(maxbrute[0])

