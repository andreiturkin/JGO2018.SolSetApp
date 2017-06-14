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
#   theta_ai_min    10              0               0
#   theta_ai_max    350             300             360
#   theta_bi_min    10              30              0
#   theta_bi_max    340             360             330
#-----------------------------------------------------------------------
#   * see Gosselin C., Jean M. - Determination of the workspace of planar
#   parallel manipulators with joint limits.
########################################################################
#
# angle = theta_ai - theta_bi + pi
#

class Phi(object):
    def __init__(self):
        # Initialize the parameters of the planar parallel robot
        # A constraints
        # self.__x_a = [-15.0, 15.0, 0.0]
        # self.__y_a = [-5.0*sqrt(3.0), -5.0*sqrt(3.0), 10.0*sqrt(3.0)]
        self.__a = []
        self.__a.append(np.array([-15.0, -5.0*sqrt(3.0)]))
        self.__a.append(np.array([15.0, -5.0*sqrt(3.0)]))
        self.__a.append(np.array([0.0, 10.0*sqrt(3.0)]))
        # B constraints
        # self.__x_b = [-5.0, 5.0, 0.0]
        # self.__y_b = [-5.0*sqrt(3.0)/3.0, -5.0*sqrt(3.0)/3.0, 10.0*sqrt(3.0)/3.0]
        self.__b = []
        self.__b.append(np.array([-5.0, -5.0*sqrt(3.0)/3.0]))
        self.__b.append(np.array([5.0, -5.0*sqrt(3.0)/3.0]))
        self.__b.append(np.array([0.0, 10.0*sqrt(3.0)/3.0]))

        self.__rho = [12.0, 27.0]
        self.__theta_a_min = [(10.0/180.0)*pi, 0.0, 0.0]
        self.__theta_a_max = [(350.0/180.0)*pi, (300.0/180.0)*pi, (360.0/180.0)*pi]
        self.__theta_b_min = [(10.0/180.0)*pi, (30.0/180.0)*pi, 0.0]
        self.__theta_b_max = [(340.0/180.0)*pi, (360.0/180.0)*pi, (330.0/180.0)*pi]

    def getRho(self):
        return self.__rho

    def getTheta_a_Bounds(self):
        return (self.__theta_a_min, self.__theta_a_max)

    def getTheta_b_Bounds(self):
        return (self.__theta_b_min, self.__theta_b_max)

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
        # Uk = [xk yk]', k in {\empty, ai, ci}
        # Uci = Uai - R(phi)
        u_c1 = self.__a[i]- self.RotL(self.__b[i], angle)
        # Rhomin^2 < Rho^2 < Rhomax^2
        # Rho^2 = norm(U - Uci)^2
        rho2i = np.linalg.norm(x-u_c1)**2
        #rho2i = (x[0]-u_c1[0])**2 + (x[1]-u_c1[1])**2
        return np.array([rho2i - self.__rho[1]**2, self.__rho[0]**2 - rho2i])

    ########################################################################
    #                               Theta a Constraints
    ########################################################################
    def getThetaA(self, x, angle, i):
        # Uk = [xk yk]', k in {\empty, ai, bi, ci, rot}
        # Ubi = U + R(phi)
        u_bi = x + self.RotL(self.__b[i], angle)

        # Rho = norm(U - Uai)
        rhoi = np.linalg.norm(u_bi-self.__a[i])

        if u_bi[1] - self.__a[i][1] >= 0.0:
            return acos((u_bi[0] - self.__a[i][0])/rhoi)
        else:
            return 2*pi - acos((u_bi[0] - self.__a[i][0])/rhoi)

    def gThetaAiBi(self, x, angle, i):
        ########################################################################
        #                               Theta a
        ########################################################################
        theta_ai = self.getThetaA(x, angle, i)

        res = []
        res.append(theta_ai - self.__theta_a_max[i])
        res.append(self.__theta_a_min[i] - theta_ai)

        ########################################################################
        #                               Theta b
        ########################################################################
        theta_bi = theta_ai - angle + pi
        if theta_bi >= 2*pi:
            theta_bi = theta_bi - 2*pi

        res.append(theta_bi - self.__theta_b_max[i])
        res.append(self.__theta_b_min[i] - theta_bi)

        return res

    def getIntl(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]
        zmin = bounds[2][0]
        zmax = bounds[2][1]
        return (interval[xmin, xmax], interval[ymin, ymax], interval[zmin, zmax])

    def getPhi(self, iBox):
        intl = self.getIntl(iBox)
        inx = [interval(i).midpoint for i in intl]
        angle = inx[2][0][0]
        x = [i[0][0] for i in inx[:2]]
        flist = ['gRhoi', 'gThetaAiBi']

        res = [r for f in flist for i in [0, 1, 2] for r in getattr(self, f)(x, angle, i)]
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

            u_c1 = self.__a[i]- Rotate(self.__b[i], xc[2])
            rho2i = NormSquared(xc[:2]-u_c1)
            ret.append(rho2i - self.__rho[1]**2)
            ret.append(self.__rho[0]**2 - rho2i)

            return ret

        def casadi_ThetaAiBi(xc, i):
            ret = []

            u_bi = xc[:2] + Rotate(self.__b[i], xc[2])
            rhoi = Norm(u_bi-self.__a[i])

            theta_ai = casadi.if_else(u_bi[1] - self.__a[i][1] >= 0.0, \
                                      casadi.arccos((u_bi[0] - self.__a[i][0])/rhoi), \
                                      2*pi-casadi.arccos((u_bi[0] - self.__a[i][0])/rhoi))
            ret.append(theta_ai - self.__theta_a_max[i])
            ret.append(self.__theta_a_min[i] - theta_ai)

            theta_bi = casadi.if_else(theta_ai - xc[2] + pi >= 2*pi, \
                                      theta_ai - xc[2] - pi, \
                                      theta_ai - xc[2] + pi)
            ret.append(theta_bi - self.__theta_b_max[i])
            ret.append(self.__theta_b_min[i] - theta_bi)

            return ret

        xc = casadi.SX.sym('xc',3)
        flist = ['casadi_Rhoi_init', 'casadi_ThetaAiBi']
        SXfs = [item for sublist in [locals()[f](xc,i) for f in flist for i in [0, 1, 2]] for item in sublist]

        maxSXfs = casadi.fmax(SXfs[17], casadi.fmax(SXfs[16],\
        casadi.fmax(SXfs[15], casadi.fmax(SXfs[14], casadi.fmax(SXfs[13], casadi.fmax(SXfs[12],\
        casadi.fmax(SXfs[11], casadi.fmax(SXfs[10], casadi.fmax(SXfs[9], casadi.fmax(SXfs[8],\
        casadi.fmax(SXfs[7], casadi.fmax(SXfs[6], casadi.fmax(SXfs[5], casadi.fmax(SXfs[4],\
        casadi.fmax(SXfs[3], casadi.fmax(SXfs[2], casadi.fmax(SXfs[1], SXfs[0])))))))))))))))))

        # Test for Taking Derivatives
        # abc = casadi.SX.sym('abc',3)
        # abcfunc = casadi.Function('abcfunc', [abc[0], abc[1], abc[2]],[abc[0]**2 + abc[1]**2 + abc[2]**2],['a','b','c'],['abcnorm2'])
        # dabcda = abcfunc.jacobian('a','abcnorm2')
        # dabcdb = abcfunc.jacobian('b','abcnorm2')
        # dabcdc = abcfunc.jacobian('c','abcnorm2')
        # norm_dabc = casadi.sqrt(dabcda(*dabcda.sx_in())[0]**2+dabcdb(*dabcdb.sx_in())[0]**2+dabcdc(*dabcdc.sx_in())[0]**2)

        phifunc = casadi.Function('phifunc', [xc[0], xc[1], xc[2]], [maxSXfs], ['x','y','z'], ['phi'])
        phi_dx = phifunc.jacobian('x','phi')
        phi_dy = phifunc.jacobian('y','phi')
        phi_dz = phifunc.jacobian('z','phi')
        dphi_norm = casadi.sqrt(phi_dx(*phi_dx.sx_in())[0]**2 + \
                                phi_dy(*phi_dy.sx_in())[0]**2 + \
                                phi_dz(*phi_dz.sx_in())[0]**2)
        dphi_norm_func = casadi.Function('dphinormfunc', [xc[0], xc[1], xc[2]], [dphi_norm], ['x','y','z'], ['normnablaphi'])
        dphi_norm_func_inv = casadi.Function('dphinorminvfunc', [xc[0], xc[1], xc[2]], [-dphi_norm], ['x','y','z'], ['normnablaphiinv'])
        dphi_norm_func_inv_jac_x = dphi_norm_func_inv.jacobian('x', 'normnablaphiinv')
        dphi_norm_func_inv_jac_y = dphi_norm_func_inv.jacobian('y', 'normnablaphiinv')
        dphi_norm_func_inv_jac_z = dphi_norm_func_inv.jacobian('z', 'normnablaphiinv')

        # Finding a max value for every norm in the list
        def iFunc(inx, *func):
            funcret = func[0].call([inx[0],inx[1],inx[2]])
            return funcret[0]

        # def iFuncInv(inx, *func):
        #     funcret = func[0].call([inx[0],inx[1],inx[2]])
        #     return -funcret[0]

        def idFuncInv(inx, *dfunc):
            ret_x = dfunc[1].call([inx[0],inx[1],inx[2]])[0]
            ret_y = dfunc[2].call([inx[0],inx[1],inx[2]])[0]
            ret_z = dfunc[3].call([inx[0],inx[1],inx[2]])[0]
            return np.array([ret_x, ret_y, ret_z])

        ret = []
        intl = self.getIntl(iBox)
        xmin, xmax = intl[0][0][0], intl[0][0][1]
        ymin, ymax = intl[1][0][0], intl[1][0][1]
        zmin, zmax = intl[2][0][0], intl[2][0][1]
        steplength = 10.0
        xstep = (xmax-xmin)/steplength
        ystep = (ymax-ymin)/steplength
        zstep = (zmax-zmin)/steplength
        rranges = (slice(xmin, xmax, xstep), slice(ymin, ymax, ystep), slice(zmin, zmax, zstep))
        maxbrute = optimize.brute(iFunc, rranges, args=[dphi_norm_func_inv], full_output=True, finish=None)
        bnds = ((xmin, xmax), (ymin, ymax), (zmin, zmax))
        ret = optimize.minimize(iFunc, maxbrute[0], args=tuple([dphi_norm_func_inv, dphi_norm_func_inv_jac_x, dphi_norm_func_inv_jac_y, dphi_norm_func_inv_jac_z]), \
                                jac=idFuncInv, method='L-BFGS-B', bounds=bnds)
        # Finding a max value for these value,
        # which is a Lipschitz constant for phi function
        return iFunc(ret.x, dphi_norm_func)

    def GlobOptExt(self, iBox):
        def getRes(inx):
            x = inx[0:2]
            angle = inx[2]
            flist = ['gRhoi', 'gThetaAiBi']
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
        zmin, zmax = intl[2][0][0], intl[2][0][1]
        steplength = 10.0
        xstep = (xmax-xmin)/steplength
        ystep = (ymax-ymin)/steplength
        zstep = (zmax-zmin)/steplength
        rranges = (slice(xmin, xmax, xstep), slice(ymin, ymax, ystep), slice(zmin, zmax, zstep))
        minbrute = optimize.brute(PhiFunc, rranges, full_output=True, finish=None)
        maxbrute = optimize.brute(PhiFuncInv, rranges, full_output=True, finish=None)

        return minbrute[1], PhiFunc(maxbrute[0]), getRes(minbrute[0]), getRes(maxbrute[0])

