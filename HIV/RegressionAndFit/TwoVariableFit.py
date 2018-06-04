import numpy as np
import scipy
from scipy.optimize import curve_fit
import csv
import math
import re

# inputs
testSetFileName = 'decideBins.txt'
pValueMatrixFileName = 'matrix.csv'
outFileName = 'out.csv'

# constants
THRESHOLD = 0.001
PATIENT = 0
DT = 1
POS_START = 2
DIGITS_PATTERN = re.compile('\d+')
TIME_THRESH = 380
ZERO_VAL_STRING = '0.001'

# compute proximal positions ahead of time
proxPosDict = {} # <int, [int pos, float dist]>

class TimePoint:
    def __init__(self, dT, patient):
        self.dT = dT
        self.patient = patient
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

def binConvertNum(num):
    if num > THRESHOLD:
        return 1
    else:
        return 0

# parse time points and binary convert volatility values
def parseTimePoints():
    with open(testSetFileName, 'r') as testSetFile:
        isFirstRow = True
        reader = csv.reader(testSetFile)
        rowCount = 0
        for row in reader:
            rowCount += 1
            if isFirstRow:
                posIndices = parsePositionIndex(row)
                allTimePoints = [] # list of time point object
                isFirstRow = False
                continue
            if len(row) == 0:
                continue
            timePoint = TimePoint(int(row[DT]), row[PATIENT])
            for key in posIndices:
                timePoint.prePos[posIndices[key]] = row[key + 1]
                timePoint.postPos[posIndices[key]] = row[key]
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
        outFile.write("Position,Temporal Persistance Weight,Spatial Persistance Weight, Total #TimePoints, #False Estimates, #Positive Mismatches, #Negative Mismatches\n")
    
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
        
        # set up fit data and do the fit
        x = []
        y = []
        z = []
        for tp in allTimePoints:
            # calculate next z
            z.append(binConvert(tp.postPos.get(posKey)))
            # calculate next y
            sigmaProx = 0.0
            for prox in allPos[posKey].proxPos:
                sigmaProx += binConvert(tp.prePos[prox]) * (1 / tp.dT) * temporalDiff(prox) * spatialDiff(adjList[prox][posKey])
            y.append(sigmaProx)
            # calculate next x
            x.append(binConvert(tp.prePos[posKey]) * (1 / tp.dT) * temporalDiff(posKey))
        
        fitResult = curve_fit(func, (x, y), zs)
        print(fitResult)
        temporalWeight = fitResult[0][0];
        spatialWeight = fitResult[0][1];
         
        # after computations are done, analyze the errors (calculate z using x and y)
        zEstimate = []
        for i in range(0, len(x)):
            zEstimate.append(func((x[i], y[i]), temporalWeight, spatialWeight))
        total = 0
        totalMismatches = 0
        positiveMismatches = 0
        negativeMismatches = 0
        for actual, estimated in zip(z, zEstimate):
            actual = binConvertNum(actual)
            estimated = binConvertNum(estimated)
            
            # do stats
            total += 1
            if actual != estimated:
                totalMismatches += 1
                if estimated > actual:
                    positiveMismatches += 1
                else:
                    negativeMismatches += 1
        
        # prepare fit results to be written to file
        nextLine = ""
        nextLine += str(posKey)
        nextLine += "," + str(temporalWeight)
        nextLine += "," + str(spatialWeight)
        
        # prepare error analysis results to outfile
        nextLine += "," + str(total)
        nextLine += "," + str(totalMismatches)
        nextLine += "," + str(positiveMismatches)
        nextLine += "," + str(negativeMismatches)
        
        # write output
        with open(outFileName, 'a') as outFile:
            outFile.write(nextLine)
        
        # write volatility values and prediction results right next to them
        with open("resultCompare.csv", 'w') as compareFile:
            writer = csv.writer(compareFile, lineterminator = '\n')
            
            # write header
            firstRow = ['patient', 'dT']
            for p in allPos[posKey].proxPos:
                firstRow.append(str(p) + "-")
                firstRow.append(str(p) + '+')
            firstRow.append('actual')
            firstRow.append('estimated')
            firstRow.append('proxmal contribution (everything except weights)')
            firstRow.append('presistence contribution (everything except weights)')
            writer.writerow(firstRow)
            
            # writer body
            z = []
            for tp in allTimePoints:
                z.append(tp.postPos.get(posKey))
                
            count = -1
            for tp in allTimePoints:
                
                # write tau values (spatial)
                tauRow = ['Tau','',]
                for p in allPos[posKey].proxPos:
                    tauRow.append(temporalDiff(p))
                    tauRow.append('')
                writer.writerow(tauRow)
                
                # write deltas
                deltaRow = ['delta','',]
                for p in allPos[posKey].proxPos:
                    deltaRow.append(spatialDiff(adjList[p][posKey]))
                    deltaRow.append('')
                writer.writerow(deltaRow)
                
                # write input data
                count += 1
                nextRow = [tp.patient, tp.dT]
                for p in allPos[posKey].proxPos:
                    nextRow.append(tp.prePos[p])
                    nextRow.append(tp.postPos[p])
                nextRow.append(z[count])
                nextRow.append(zEstimate[count])
                nextRow.append(y[count])
                nextRow.append(x[count])
                writer.writerow(nextRow)
                writer.writerow([])
                    