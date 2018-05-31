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

# compute proximal positions ahead of time
proxPosDict = {} # <int, [int pos, float dist]>

class timePoints:
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
        pos = int(digitsPattern.match(row[i]).group(0))
        indices[i] = pos # index -> vn
    return indices

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
            timePoint = TimePoint(int(row[DT]))
            for key in posIndices:
                timePoint.prePos[posIndices[key]] = row[key + 1]
                timePoint.prePos[posIndices[key]] = row[key]
    return allTimePoints

# read a distance matrix file, 
def readDistMatrix():
    with open(pValueMatrixFileName, 'r') as matrixFile:
        isFirstRow = True
        reader = csv.reader(matrixFile)
        indToPosDict = {}
        adjList = {}
        colInd = 0
        for line in reader:
            if isFirstRow:
                for i in range(1, len(line)):
                    indToPosDict[i] = int(line[i])
                isFirstRow = False
                for key in posToIndDict: # initialize adjacency list
                    adjList[key] = {}
                continue
            for i in range(1, len(line)):
                adjList[indToPosDict[colInd]][indToPosDict[i]] = float[row[i]]
                adjList[indToPosDict[i]][indToPosDict[colInd]] = float[row[i]]
            colInd += 1
        return [adjList, indToPosDict]

# diff functions
def spatialDiff(pValue):
    return 2 / (1 + math.exp(-(pValue  - 2.5) * 0.3)) - 1

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
    elif pos == 330:
        return 0.03
    elif pos == 415:
        return 0.03
    elif pos == 495:
        return 0
    elif pos == 519:
        return 0.003
    elif pos == 770:
        return 0.01
    else:
        raise Exception("not yet implemented")

# run the two variable fit
if __name__ == "__main__":
    
    # clear output file and write headers
    with open(outFileName, 'w') as outFile:
        outFile.write("Position,Temporal Persistance Weight,Spatial Persistance Weight\n")
    
    # parse adjList for pValue distance
    [adjList, indToPosDict] = readDistMatrix()
    
    # parse all time points
    allTimePoints = parseTimePoints()
    
    # set up proximal positions based on adjList
    allPos = {}
    for key in adjList:
        if key == 332: #TODO remove this condition later, for now we doing 332 only
            allPos[key] = Position(key)
    for key in allPos:
        for proxKey in adjList[key]:
            if adjList[key][adjKey] < THRESHOLD:
                allPos[key].proxPos.append(adjKey)
    
    # compute independent variables for each position
    for pos in allPos:
        x = []
        y = []
        z = []        
        for tp in allTimePoints:
            for posKey in allPos:
                # calculate next z
                z.append(tp.postPos.get(posKey))
                # calculate next y
                sigmaProx = 0.0
                for prox in allPos[posKey].proxPos:
                    sigmaProx += tp.prePos[prox] * (1 / tp.dT) * temporalDiff(prox) * (1 / spatialDiff(adjList[prox][posKey]))
                y.append(sigmaProx)
                # calculate next x
                x.append(tp.prePos[posKey] * (1 / tp.dT) * temporalDiff(posKey))
                
                # initial guesses for a, b:
                # p0 = 1, 2
                
                # compute fit parameters and write to file output   
                fitResult = curve_fit(func, (x,y), z, p0)
                with open(outFileName, 'a') as outFile:
                    nextLine = ""
                    nextLine += str(pos.posNum)
                    nextLine += "," + str(fixResult[0][0])
                    nextLine += "," + str(fixResult[0][1]) + "\n"
                    outFile.write(nextLine)