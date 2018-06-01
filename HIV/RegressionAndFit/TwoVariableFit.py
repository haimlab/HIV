import numpy as np
import scipy
from scipy.optimize import curve_fit
import csv
import math
import re

# inputs
testSetFileName = 'testSet.csv'
pValueMatrixFileName = 'matrix.csv'
outFileName = 'out.csv'

# constants
THRESHOLD = 0.001
DT = 1
POS_START = 2
DIGITS_PATTERN = re.compile('\d+')
TIME_THRESH = 380
ZERO_VAL_STRING = '0.001'

# compute proximal positions ahead of time
proxPosDict = {} # <int, [int pos, float dist]>

class TimePoint:
    def __init__(self, dT):
        self.dT = dT
        self.prePos = {}
        self.postPos = {} # <int, float>
        
class Position:
    def __init__(self, posNum):
        self.posNum = posNum
        self.timePoints = {} # dt -> [vn -1, vn]
        self.proxPos = []
    def addTimePoint(self, dT, Vpost, Vpre):
        self.timePoints.append(toAdd)

# the function to run the fit on
def func(X, a, b):
    x, y = X
    return a * x + b * y

# get the index of the positions, assumes the first row as a list
def parsePositionIndex(row):
    indices = {}
    for i in range(POS_START, len(row), 2):
        if row[i].strip() != '':
            pos = int(DIGITS_PATTERN.match(row[i]).group(0))
            indices[i] = pos # index -> vn
    return indices

# if the string represent a number, convert it to 0 if it equals 0.001 and 1
# otherwise, raise an exception if myString is not float
def binConvert(myString):
    if not myString.replace('.', '').replace('E', '').replace('-', '').isdigit():
        raise Exception('The string is not float format')
    if myString == ZERO_VAL_STRING:
        return 0
    else:
        return 1


# parse time points and binary convert volatility values
def parseTimePoints():
    with open(testSetFileName, 'r') as testSetFile:
        isFirstRow = True
        reader = csv.reader(testSetFile)
        for row in reader:
            if isFirstRow:
                posIndices = parsePositionIndex(row)
                allTimePoints = [] # list of time point object
                isFirstRow = False
                continue
            if len(row) == 0:
                continue
            timePoint = TimePoint(int(row[DT]))
            for key in posIndices:
                timePoint.prePos[posIndices[key]] = binConvert(row[key + 1])
                timePoint.postPos[posIndices[key]] = binConvert(row[key])
            if timePoint.dT < TIME_THRESH:
                allTimePoints.append(timePoint)
    return allTimePoints

# read a distance matrix file, 
def readDistMatrix():
    with open(pValueMatrixFileName, 'r') as matrixFile:
        isFirstRow = True
        reader = csv.reader(matrixFile)
        indToPosDict = {}
        adjList = {}
        colInd = 1
        for line in reader:
            if isFirstRow:
                for i in range(1, len(line)):
                    indToPosDict[i] = int(line[i])
                isFirstRow = False
                for key in indToPosDict: # initialize adjacency list
                    adjList[indToPosDict[key]] = {}
                continue
            for i in range(1, len(line)):
                if line[i].strip() != '':
                    adjList[indToPosDict[colInd]][indToPosDict[i]] = float(line[i])
                    adjList[indToPosDict[i]][indToPosDict[colInd]] = float(line[i])
            colInd += 1
        return adjList, indToPosDict

# diff functions
def spatialDiff(pValue):
    if pValue == 0:
        pValue = 10e-9
    pValue = -math.log(pValue, 10)
    result = 2 / (1 + math.exp(-(pValue  - 2.5) * 0.3)) - 1
    if result < 0:
        result = 0
    return result

# return the temporal coefficient of a specific 
def temporalDiff(pos):
    if pos == 117:
        return 0
    elif pos == 135:
        return 0.05
    elif pos == 147:
        return 0.05
    elif pos == 175:
        return 0
    elif pos == 181:
        return 0.01
    elif pos == 295:
        return 0.03
    elif pos == 299:
        return 0
    elif pos == 318:
        return 0
    elif pos == 327:
        return 0
    elif pos == 332:
        return 0.35
    elif pos == 330 or pos == 334:
        return 0.03
    elif pos == 415:
        return 0.03
    elif pos == 495:
        return 0
    elif pos == 519:
        return 0.003
    elif pos == 779:
        return 0.01
    else:
        raise Exception("not yet implemented")

# run the two variable fit
if __name__ == "__main__":
    
    # clear output file and write headers
    with open(outFileName, 'w') as outFile:
        outFile.write("Position,Temporal Persistance Weight,Spatial Persistance Weight\n")
    
    # parse adjList for pValue distance
    adjList, indToPosDict = readDistMatrix()
    
    # parse all time points
    allTimePoints = parseTimePoints()
    
    # set up proximal positions based on adjList
    allPos = {}
    for key in adjList:
        if key == 332: #TODO remove this condition later, for now we doing 332 only
            allPos[key] = Position(key)
    for key in allPos:
        for proxKey in adjList[key]:
            if adjList[key][proxKey] < THRESHOLD:
                allPos[key].proxPos.append(proxKey)
    
    # compute independent variables for each position
    for posKey in allPos:
        x = []
        y = []
        z = []
        for tp in allTimePoints:
            # calculate next z
            z.append(tp.postPos.get(posKey))
            # calculate next y
            sigmaProx = 0.0
            for prox in allPos[posKey].proxPos:
                sigmaProx += tp.prePos[prox] * (1 / tp.dT) * temporalDiff(prox) * spatialDiff(adjList[prox][posKey])
            y.append(sigmaProx)
            # calculate next x
            x.append(tp.prePos[posKey] * (1 / tp.dT) * temporalDiff(posKey))
            
        # initial guesses for a, b:
        # p0 = 1, 2
            
        # compute fit parameters and write to file output   
        print(len(x))
        print(len(y))
        print(len(z))
        
        fitResult = curve_fit(func, (x,y), z)
        with open(outFileName, 'a') as outFile:
            nextLine = ""
            nextLine += str(posKey)
            nextLine += "," + str(fitResult[0][0])
            nextLine += "," + str(fitResult[0][1]) + "\n"
            outFile.write(nextLine)