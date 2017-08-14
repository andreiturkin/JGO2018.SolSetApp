import os, glob, gc
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
from JGO2017rev_Example2 import Example2_Lipz
from JGO2017rev_Example2 import Example2_LipzHcs
from JGO2017rev_Example3 import Example3
from JGO2017rev_NUC import Box

def CreateOrClear(DirName, bClean):
    if not os.path.isdir(DirName):
        os.makedirs(DirName)
    else:
        if bClean:
            for f in glob.glob(DirName + '*.*'):
                os.remove(f)

def JGOExample1(iDelta, ShowRes=False, Zoom=False, ZoomBox=None):

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
            if Zoom and ZoomBox:
                System.SaveResults('./Images/{}_Delta_{}_Zoom.pdf'.format(flist[i], iDelta), False, Zoom, ZoomBox)
            else:
                System.SaveResults('./Images/{}_Delta_{}.pdf'.format(flist[i], iDelta), AddRings=True)

        if exectime:
            retdata.append(['{}'.format(iDelta), '{}'.format(System.getResIterations()), '{}'.format(round(exectime,4))])
        else:
            retdata.append(['{}'.format(iDelta), 'Unknown', 'Unknown'])
    return [item for dt in retdata for item in dt]

def JGOExample2(iDelta, ShowRes=False, Zoom=False, ZoomBox=None):
    # Define a system of inequalities
    flist = ['Example2_Lipz', 'Example2_LipzHcs']
    retdata = []
    for i, funcname in enumerate(flist):
        System = globals()[flist[i]](iDelta, ShowCovPrc=False)
        fName = './Dumps/{}_{}.p'.format(funcname, iDelta)
        exectime = None
        if System.isFileExist(fName):
            System.LoadSolution(fName)
        else:
            maxLevels = 64
            t = timeit.Timer(lambda: System.getSolution(maxLevels))
            exectime = t.timeit(number=1)
            System.SaveSolution(fName)

        if ShowRes:
            if Zoom and ZoomBox:
                System.SaveResults('./Images/{}_Delta_{}_Zoom.pdf'.format(funcname, iDelta), Zoom, ZoomBox)
            else:
                System.SaveResults('./Images/{}_Delta_{}.pdf'.format(funcname, iDelta))
        retdata.append(['{}'.format(iDelta), '{}'.format(System.getResIterations()), '{}'.format(round(exectime,4))] \
                       if exectime else ['{}'.format(iDelta), 'Unknown', 'Unknown'])
    return [item for dt in retdata for item in dt]

def JGOExample3(iDelta, iAngle, ShowRes=False):
    # Define a system of inequalities
    System = Example3(idelta=iDelta, iangle=iAngle, ShowCovPrc=False)
    fName = './Dumps/Example3_{}_{}.p'.format(iDelta, round(iAngle*180.0/pi))
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
        System.SaveResults('./Images/Example3_Delta_{}_Angle_{}.pdf'.format(iDelta, round(iAngle*180.0/pi)), AddRings=True)
    return ['{}'.format(iDelta), '{}'.format(System.getResIterations()), '{}'.format(round(exectime,4))] \
           if exectime else ['{}'.format(iDelta), 'Unknown', 'Unknown']

if __name__ == '__main__':
    ImageDir = 'Images/'
    DumpDir = 'Dumps/'
    # Cleanning up the directories
    CreateOrClear(ImageDir, True)
    CreateOrClear(DumpDir, True)

    print'###############################################################################'
    print'             Approximating The Solution Set of Nonlinear inequalities            '
    print'###############################################################################'
    print'                            Example 1. Table 1'
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
    print'                            Example 1. Figure 1'
    JGOExample1(0.06, True)

    # Uncomment the following 6 lines if you would to get the result as it is in the Figure 2
    # print'                            Example 1. Figure 2'
    # print'                           (it will take a while)'
    # gc.collect()
    # JGOExample1(0.0001, True)
    # gc.collect()
    # JGOExample1(0.0001, True, True, Box((1.0-0.002, -0.002), (0.004, 0.004)))

    print'###############################################################################'
    print'                            Example 2. Table 2'
    Example2_deltas = [0.25, 0.1, 0.05, 0.03, 0.018, 0.015, 0.01, 0.005]
    pbar = progressbar.ProgressBar(maxval=len(Example2_deltas), \
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    pbar.start()
    rows = []
    rows.append([' ','Analytically calculated L',' ',' ','Numerically calculated L',' '])
    rows.append(['Delta', 'Iterations', 'Time', 'Delta', 'Iterations', 'Time'])
    for i, d in enumerate(Example2_deltas):
        rows.append(JGOExample2(d, False))
        pbar.update(i+1)
    pbar.finish()
    table_instance = SingleTable(rows)
    table_instance.inner_heading_row_border = False
    print(table_instance.table)
    print'                             Example 2. Figure 4'
    JGOExample2(0.01, True)

    print'###############################################################################'
    print'                            Example 3. Table 4'
    Example3_deltas = [0.1, 0.07, 0.05, 0.04, 0.025, 0.02, 0.012, 0.009]
    bar = progressbar.ProgressBar(maxval=len(Example3_deltas), \
          widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    rows = []
    rows.append([' ','Locally Calculated L',' '])
    rows.append(['Delta', 'Iterations', 'Time'])
    for i, d in enumerate(Example3_deltas):
        rows.append(JGOExample3(iDelta=d, iAngle=(50.0/180.0)*pi, ShowRes=False))
        bar.update(i+1)
    bar.finish()
    table_instance = SingleTable(rows)
    table_instance.inner_heading_row_border = False
    print(table_instance.table)

    print'                             Example 3. Figure 5'
    angles = [(a/180.0)*pi for a in [50, 80, 120, 140]]
    bar = progressbar.ProgressBar(maxval=len(angles), \
          widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    for i, angle in enumerate(angles):
        JGOExample3(iDelta=0.1, iAngle = angle, ShowRes=True)
        bar.update(i+1)
    bar.finish()

    print'###############################################################################'
    print' Done!'
    print' Images are in {}'.format(ImageDir)
    print' All coverings were saved as Pickle files to {} folder'.format(DumpDir)
    print'###############################################################################'
