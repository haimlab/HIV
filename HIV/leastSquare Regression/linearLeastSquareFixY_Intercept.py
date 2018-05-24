import csv
import numpy
import math

# TODO stream line application of rules

# inputs
binAnalysisOutFileName = 'C:\\Users\\rdong6\\Desktop\\least square dev\\out consec.csv'
decideBinsFileName = 'C:\\Users\\rdong6\\Desktop\\least square dev\\decideBins.csv'
outFileName = 'C:\\Users\\rdong6\\Desktop\\least square dev\\output.csv'
binIntervals = [(1, 61), (62, 200), (201, 670), (671, 4433)] # in order
threshold = 0.05

# constants
dt = 1
pValStart = 1
nonSignificant = 1

# apply rule 1 to fit data object
# rule 1: if p-value of the first bin is not significant, then every bin in that
# position becomes non-significant
def applyRule1(fitData):
    if fitData.yData[0] > threshold:
        length = len(fitData.yData)
        for i in range(0, length):
            fitData.yData[i] = nonSignificant

# apply rule 2 to fit data object
# rule 2: if there exists two consecutive p-Values that are both insignificant,
# all following p-Values also become insignificant
def applyRule2(fitData):
    length = len(fitData.yData)
    if length < 2:
        raise Exception('need at least 2 data points for rule 2')
    pre = fitData.yData[0] > threshold
    for i in range(1, length):
        cur = fitData.yData[i] > threshold
        if pre and cur:
            for j in range(i - 1, length):
                fitData.yData[j] = 1
            break
        pre = cur

# takes a list of x coordinates, and shift them such that the first point lies 
# exactly on the origin
def normalizeToOrigin(myList):
    if len(myList) == 0:
        raise Exception('list is empty')
    myList.sort()
    shiftDist = myList[0]
    for i in range(0, len(myList)):
        myList[i] -= shiftDist
    
# perfrom linear least square fit with fixed intercept
# X: x coordinates of all sample points
# Y: y coordinates of all sample points
# intcpt: the y-intercept to fix to
def linearFixIntcptLSF(X, Y, intcpt):
    guardDim(X, Y)
    xSum = 0
    xSquaredSum = 0
    xYSum = 0
    for i in range(0, len(X)):
        xSum += X[i]
        xSquaredSum += X[i] ** 2
        xYSum += X[i] * Y[i]
    return (xYSum - intcpt * xSum) / xSquaredSum

# guard against dimension mismatch
def guardDim(X, Y):
    if len(X) != len(Y):
        raise Exception('dimension mismatch')

class FitData:
    def __init__(self, medians):
        self.xData = medians
        self.yData = []
    def addData(self, toAdd):
        self.yData.append(toAdd)
    def calcFit(self):
        self.slope = linearFixIntcptLSF(self.xData, self.yData, self.yData[0])
        self.yIntcpt = self.yData[0]
        try:
            self.xIntcpt = (-math.log(threshold, 10) - self.yIntcpt) / self.slope
        except ZeroDivisionError:
            self.xIntcpt = 'slope = 0'
    def getFitEquation(self):
        return 'y = ' + str(self.slope) + 'x' + ' + ' + str(self.yIntcpt)

class Bin:
    def __init__(self, interval):
        if len(interval) != 2:
            raise Exception('incorrect format for interval')
        self.start = int(interval[0])
        self.end = int(interval[1])
        self.timePoints = []
    def containsDT(self, dt):
        dt = int(dt)
        return dt <= self.end and dt >= self.start
    def addDt(self, dt):
        dt = int(dt)
        if self.containsDT(dt):
            self.timePoints.append(dt)
    def median(self):
        length = len(self.timePoints)
        if length == 0:
            raise Exception('cannot compute median of empty bin')
        self.timePoints.sort()
        if length % 2 == 0: # even length array
            a = length // 2
            b = length // 2 - 1
            return (self.timePoints[a] + self.timePoints[b]) / 2
        else: # odd length array
            return self.timePoints[length // 2]            

# do linear fit for output file of bin-analysis
# compute the medians
dtDict = {}
bins = []
for interval in binIntervals:
    bins.append(Bin(interval))
with open(decideBinsFileName, 'r') as binInFile:
    reader = csv.reader(binInFile)
    firstRow = True
    for row in reader:
        if firstRow:
            firstRow = False
            continue
        t = int(row[dt])
        for b in bins:
            b.addDt(t)
medians = []
for b in bins:
    medians.append(b.median())
normalizeToOrigin(medians)

# read in data for least square fit for all positions
allPos = []
allPosFitData = {}
with open(binAnalysisOutFileName, 'r') as binOutFile:
    reader = csv.reader(binOutFile)
    firstRow = True
    for row in reader:
        if firstRow: # get all positions
            row = row[1:]
            for p in row:
                allPos.append(int(p))
                allPosFitData[int(p)] = FitData(medians)
            firstRow = False
            continue
        row = row[pValStart:]
        for i in range(0, len(allPos)):
            allPosFitData[allPos[i]].addData(float(row[i]))

# apply rule 1 and 2 to all fit data
for pos in allPos:
    applyRule1(allPosFitData[pos])
    applyRule2(allPosFitData[pos])

# correct p values w.r.t 1 using -log_10(p)
for key in allPosFitData:
    yData = allPosFitData[key].yData
    length = len(yData)
    for i in range(0, length):
        yData[i] = -math.log(yData[i], 10)

# compute fit parameters for each position
for key in allPosFitData:
    allPosFitData[key].calcFit()

# write output
with open(outFileName, 'w') as outFile, \
     open(binAnalysisOutFileName, 'r') as oldOutFile:
    # put a copy of input files contents for convinience
    contents = oldOutFile.read()
    contents = contents.strip()
    outFile.write(contents + "\n")
    # write the medians and fit data after applying rules for convinience
    outFile.write('bin medians, fit data after applying rules + \n')
    numBins = len(medians)
    for i in range(0, numBins):
        outFile.write(str(medians[i]) + ',')
        for pos in allPos:
            outFile.write(str(allPosFitData[pos].yData[i]) + ',')
        outFile.write('\n')
    # write the fit results
    outFile.write('fit equation,')
    for key in allPosFitData:
        outFile.write(allPosFitData[key].getFitEquation() + ',')
    outFile.write('\n')
    outFile.write('loose significance at,')
    for key in allPosFitData:
        outFile.write(str(allPosFitData[key].xIntcpt) + ',')