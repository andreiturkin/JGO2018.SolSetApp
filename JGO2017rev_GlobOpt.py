from JGO2017rev_NUC import Box

class GlobalOptimization(object):
    def __init__(self, iFunc, igetLipValFunc, ieps=0):
        # The Initial Parameters
        self.__eps = ieps
        # These two functions are used to get a minorant and majorant
        # for the function is used
        self.__Func = iFunc
        self.__getLipValFunc = igetLipValFunc

    ############################################################################################
    # Global Optimization
    #
    ############################################################################################
    def getCurrentMinVal(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        return self.__Func(((xmin+xmax)/2.0, (ymin+ymax)/2.0))-(self.__getLipValFunc(iBox)/2.0)*iBox.getDiam()

    def getCurrentMaxVal(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        return self.__Func(((xmin+xmax)/2.0, (ymin+ymax)/2.0))+(self.__getLipValFunc(iBox)/2.0)*iBox.getDiam()

    def getCurrentCntrVal(self, iBox):
        bounds = iBox.getBounds()
        xmin = bounds[0][0]
        xmax = bounds[0][1]
        ymin = bounds[1][0]
        ymax = bounds[1][1]

        return self.__Func(((xmin+xmax)/2.0, (ymin+ymax)/2.0))

    ############################################################################################
    # Global Optimization
    # with Predefined Accuracy (eps)
    ############################################################################################
    def getMinVal(self, iBox):
        curBoxes = [iBox]
        lCntrs = [self.getCurrentCntrVal(iBox)]
        tempBoxes = []
        while curBoxes:
            for box in curBoxes:
                minorant = self.getCurrentMinVal(box)
                if minorant >= min(lCntrs) - self.__eps:
                    continue

                lCntrs.append(self.getCurrentCntrVal(box))

                box1, box2 = box.Split()
                tempBoxes.append(box1)
                tempBoxes.append(box2)

            curBoxes = tempBoxes
            tempBoxes = []

        minVal = min(lCntrs)
        lCntrs = []

        return minVal

    def getMaxVal(self, iBox):
        curBoxes = [iBox]
        lCntrs = [self.getCurrentCntrVal(iBox)]
        tempBoxes = []
        while curBoxes:
            for box in curBoxes:
                majorant = self.getCurrentMaxVal(box)
                if majorant <= max(lCntrs) + self.__eps:
                    continue

                lCntrs.append(self.getCurrentCntrVal(box))

                box1, box2 = box.Split()
                tempBoxes.append(box1)
                tempBoxes.append(box2)

            curBoxes = tempBoxes
            tempBoxes = []

        maxVal = max(lCntrs)
        lCntrs = []

        return maxVal