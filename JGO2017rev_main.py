import timeit
import datetime
import matplotlib
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
matplotlib.use('TkAgg')

import numpy as np

# A Parallel Robot
from JGO2017rev_Example2 import Example2

def GetWorkspace2(iDelta, ShowRes=False):
    # Define a system of inequalities
    bIntl = False
    System = Example2(iDelta, ShowCovPrc=False)
    fName = './Dumps/Example2_Tree_{}_lipz.p'.format(iDelta)
    print 'The results for {}:'.format(iDelta)
    if System.isFileExist(fName):
        print 'A file with a precalculated solution was found: {}'.format(fName)
        System.LoadSolution(fName)
    else:
        maxLevels = 64
        t = timeit.Timer(lambda: System.getSolution(maxLevels))
        exectime = t.timeit(number=1)
        print 'Execution time: {}\n'.format(exectime)
        print '...Saving Data to {}'.format(fName)
        System.SaveSolution(fName)

    if ShowRes:
        System.SaveResults('./Images/Example2_{0}__{1:02d}_{2:02d}_{3:02d}_covering_{4}.pdf'.format(datetime.date.today(), \
                                                           datetime.datetime.now().hour,\
                                                           datetime.datetime.now().minute,\
                                                           datetime.datetime.now().second,\
                                                           iDelta), AddRings=False)

if __name__ == '__main__':
    print'\n#############################################################################'
    print'                            Parallel Robot                                     '
    print'#############################################################################\n'
    #Use different delta values to get the workspace
    deltas = [0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.06, 0.035, 0.03, 0.018, 0.01, 0.009]

    # Get the workspace
    GetWorkspace2(0.1, True)

    print'\n#############################################################################'
    print'                                 Done!                                         '
    print'#############################################################################\n'
