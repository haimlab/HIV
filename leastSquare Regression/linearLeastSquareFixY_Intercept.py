import csv
import numpy

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
# inputs
binAnalysisOutFileName = ''
decideBinsFileName = 'decideBins.txt'
outFileName = ''
binIntervals = [(1, 61), (62, 200), (201, 670), (671, 4433)] # in order

# constants
dt = 1
pValStart = 1

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
    
