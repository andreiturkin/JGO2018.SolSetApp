import os, glob
import timeit
import datetime
import progressbar
from terminaltables import SingleTable
import matplotlib
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
matplotlib.use('TkAgg')
import warnings
warnings.filterwarnings("ignore")

import numpy as np
from numpy import pi

# A Parallel Robot
from JGO2017rev_Example1 import Example1_GlobOpt
from JGO2017rev_Example1 import Example1_AppxGlobL
from JGO2017rev_Example1 import Example1_AppxLocL
from JGO2017rev_Example2 import Example2
from JGO2017rev_Example3 import Example3

def CreateOrClear(DirName):
    if not os.path.isdir(DirName):
        os.makedirs(DirName)
    else:
        for f in glob.glob(DirName + '*.*'):
            os.remove(f)

def JGOExample1(iDelta, ShowRes=False):

    flist = ['Example1_GlobOpt', 'Example1_AppxGlobL', 'Example1_AppxLocL']
    retdata = []
    for i in range(3):
        System = globals()[flist[i]](iDelta, ShowCovPrc=False)
        fName = './Dumps/{}_{}.p'.format(flist[i], iDelta)
        exectime = None
        if System.isFileExist(fName):
            System.LoadSolution(fName)
        else:
            maxLevels = 64
            t = timeit.Timer(lambda: System.getSolution(maxLevels))
            exectime = t.timeit(number=1)
            System.SaveSolution(fName)

        if ShowRes:
            System.SaveResults('./Images/{}_Delta_{}.pdf'.format(flist[i], iDelta), AddRings=True)
        if exectime:
            retdata.append(['{}'.format(iDelta), '{}'.format(System.getResIterations()), '{}'.format(round(exectime,4))])
        else:
            retdata.append(['{}'.format(iDelta), 'Unknown', 'Unknown'])
    return [item for dt in retdata for item in dt]

def JGOExample2(iDelta, iAngle, ShowRes=False):
    # Define a system of inequalities
    System = Example2(idelta=iDelta, iangle=iAngle, ShowCovPrc=False)
    fName = './Dumps/Example2_{}_{}.p'.format(iDelta, round(iAngle*180.0/pi))
    # print 'The results for {}:'.format(iDelta)
    exectime = None
    if System.isFileExist(fName):
        System.LoadSolution(fName)
    else:
        maxLevels = 64
        t = timeit.Timer(lambda: System.getSolution(maxLevels))
        exectime = t.timeit(number=1)
        System.SaveSolution(fName)

    if ShowRes:
        System.SaveResults('./Images/Example2_Delta_{}_Angle_{}.pdf'.format(iDelta, round(iAngle*180.0/pi)), AddRings=True)
    return ['{}'.format(iDelta), '{}'.format(System.getResIterations()), '{}'.format(round(exectime,4))] \
           if exectime else ['{}'.format(iDelta), 'Unknown', 'Unknown']

def JGOExample3(iDelta, ShowRes=False):
    # Define a system of inequalities
    System = Example3(iDelta, ShowCovPrc=False)
    fName = './Dumps/Example3_{}.p'.format(iDelta)
    exectime = None
    if System.isFileExist(fName):
        System.LoadSolution(fName)
    else:
        maxLevels = 64
        t = timeit.Timer(lambda: System.getSolution(maxLevels))
        exectime = t.timeit(number=1)
        System.SaveSolution(fName)

    if ShowRes:
        System.SaveResults('./Images/Example3_Delta_{}.pdf'.format(iDelta), AddRings=False)
    return ['{}'.format(iDelta), '{}'.format(System.getResIterations()), '{}'.format(round(exectime,4))] \
           if exectime else ['{}'.format(iDelta), 'Unknown', 'Unknown']

if __name__ == '__main__':
    ImageDir = 'Images/'
    DumpDir = 'Dumps/'
    # Cleanning up the directories
    CreateOrClear(ImageDir)
    CreateOrClear(DumpDir)

    print'###############################################################################'
    print'                            Nonuniform covering                                '
    print'###############################################################################'
    print'                             Example 1. Table 1'
    Example1_deltas = [0.08, 0.06, 0.035, 0.03, 0.018, 0.015, 0.01, 0.001]
    bar = progressbar.ProgressBar(maxval=len(Example1_deltas), \
          widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    rows = []
    rows.append([' ','Global Optimization',' ', ' ','Global L',' ', ' ','Local L',' '])
    rows.append(['Delta', 'Iterations', 'Time', 'Delta', 'Iterations', 'Time', 'Delta', 'Iterations', 'Time'])
    for i, d in enumerate(Example1_deltas):
        rows.append(JGOExample1(d, False))
        bar.update(i+1)
    bar.finish()
    table_instance = SingleTable(rows)
    table_instance.inner_heading_row_border = False
    print(table_instance.table)

    print'###############################################################################'
    print'                               Example 1. Figure 2'
    JGOExample1(0.06, True)
    print'###############################################################################'
    print'                               Example 2. Table 2'
    Example2_deltas = [0.1, 0.07, 0.05, 0.04, 0.025, 0.02, 0.012, 0.009]
    bar = progressbar.ProgressBar(maxval=len(Example2_deltas), \
          widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    rows = []
    rows.append([' ','Local L',' '])
    rows.append(['Delta', 'Iterations', 'Time'])
    for i, d in enumerate(Example2_deltas):
        rows.append(JGOExample2(iDelta=d, iAngle=(50.0/180.0)*pi, ShowRes=False))
        bar.update(i+1)
    bar.finish()
    table_instance = SingleTable(rows)
    table_instance.inner_heading_row_border = False
    print(table_instance.table)

    print'                               Example 2. Figure 5'
    angles = [(a/180.0)*pi for a in [50, 80, 120, 140]]
    bar = progressbar.ProgressBar(maxval=len(angles), \
          widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    for i, angle in enumerate(angles):
        JGOExample2(iDelta=0.1, iAngle = angle, ShowRes=True)
        bar.update(i+1)
    bar.finish()

    print'###############################################################################'
    print'                               Example 3. Figure'
    JGOExample3(0.01, True)
    print'###############################################################################'
    print' Done!'
    print' Images are in {}'.format(ImageDir)
    print' All coverings were saved as Pickle files to {} folder'.format(DumpDir)
    print'###############################################################################'
