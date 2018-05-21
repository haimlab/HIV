import csv
import numpy

# inputs
binAnalysisOutFileName = 'C:\\Users\\rdong6\\Desktop\\leastSquareTest\\out consec.csv'
decideBinsFileName = 'C:\\Users\\rdong6\\Desktop\\leastSquareTest\\decideBins.csv'
outFileName = 'C:\\Users\\rdong6\\Desktop\\leastSquareTest\\output.csv'
binIntervals = [(1, 61), (62, 200), (201, 670), (671, 4433)] # in order
threshold = 0.05

# constants
dt = 1
pValStart = 1

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

class Position:
    def __init__(self, pos):
        self.position = int(pos)
        self.pValues = []
    def addPValue(self, toAdd):
        self.pValues.append(int(toAdd))

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
            self.xIntcpt = (threshold - self.yIntcpt) / self.slope
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
            print("bin range " + str(self.start) + ' - ' + str(self.end), flush = True)
            print("appendin " + str(dt))
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
print('medians are:')
for b in bins:
    print(str(b.median()))
    medians.append(b.median())

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

# compute fit parameters for each position
for key in allPosFitData:
    allPosFitData[key].calcFit()

# write output
with open(outFileName, 'w') as outFile, \
     open(binAnalysisOutFileName, 'r') as oldOutFile:
    contents = oldOutFile.read()
    contents = contents.strip()
    outFile.write(contents + "\n")
    outFile.write('fit equation,')
    for key in allPosFitData:
        outFile.write(allPosFitData[key].getFitEquation() + ',')
    outFile.write('\n')
    outFile.write('loose significance at,')
    for key in allPosFitData:
        outFile.write(str(allPosFitData[key].xIntcpt) + ',')
    

        
        
