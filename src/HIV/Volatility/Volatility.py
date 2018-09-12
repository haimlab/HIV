# currently we ignore the time point after patient name and use years as the unit of time points

import os
from look_up_tables import lookup
import numpy as np

# all input variables
headerOffset = 4  # number head columns
site = '2G12'  # 2G12/2F5/all
isFiveD = False
lookupSource = 'bindingefficiency'
geneticDistanceFolderName = "C:\\Users\\rdong6\\Desktop\\C_5D_And_Regular_Volatility\\Genetic_Distance\\Hyphy_Output\\"
pngsTableFileName = "C:\\Users\\rdong6\\Desktop\\C_5D_And_Regular_Volatility\\Pngs_Conversion\\Outputs\\gaps_removed" \
                    "\\all_together_formatted.txt "
# end of all input variables

# other constants
distanceFileNames = os.listdir(geneticDistanceFolderName)
skipLessThanTwoVolatilities = True


# takes in two vectors, and multiply number for each dimension
def weightedEuclidean(vec1, vec2, weightVector):
    diffVector = [v1 - v2 for v1, v2 in zip(vec1, vec2)]
    if len(diffVector) != len(weightVector):
        raise Exception("Dimensions of input vectors mismatch")
    mySum = 0
    for i in range(0, len(diffVector)):
        mySum += diffVector[i] ** 2 * weightVector[i]
    return mySum ** .5


# read genetic distances as dictionary mapping a tuple of id to distance, group by time points
def readGeneticDistance(distFileName):
    genDistMatrix = {}
    with open(distFileName, 'r') as distFile:
        isFirstRow = True
        for row in distFile:
            if isFirstRow:  # skip first row
                isFirstRow = False
                continue
            row = row.split()
            genDistMatrix[(row[0], row[1])] = float(row[2])
    return genDistMatrix


# read in pngs file as dictionary mapping from id pair to Amino Acid Sequences
def readPngs(pngsFileName):
    # allPosDict maps positions to index in pngs table relative to first letter
    allPosDict, pngsTable = dict(), dict()
    with open(pngsFileName, 'r') as pngsFile:
        isFirstLine = True
        for line in pngsFile:
            line = line.strip().split()
            if isFirstLine:  # create position to index dictionary using header row
                isFirstLine = False
                ind = 0
                for pos in line[headerOffset:]:
                    allPosDict[int(pos)] = ind
                    ind += 1
                continue
            pngsTable['.'.join(line[:headerOffset])] = line[headerOffset:]
    return pngsTable, allPosDict


# takes a vaccineSite and a list of positions ordered the same way as that in pngs file. Returns the positions we want
# to calculate volatility for and their w.r.t row (including header columns)
def getNeededPositions(site, allPosDict):
    if site == '2G12':
        needed = {295, 332, 339, 392, 448}
    elif site == '2F5':
        needed = {662, 663, 664, 665, 667}
    elif site == 'all':
        needed = set([pos for pos in allPosDict])  # all present positions
    else:
        raise Exception("unidentified position set")
    neededPosInds = dict()
    for pos in allPosDict:
        if pos in needed:
            neededPosInds[pos] = allPosDict[pos]
    return neededPosInds


# calculate volatility
def calcVolatility(genDistMatrix, pngsTable, neededPosInds, isFiveD, lookupSource):
    volatilities = []
    if isFiveD:
        total = 0
        for id1, id2 in genDistMatrix:
            vec1, vec2 = [], []
            for pos in neededPosInds:
                vec1.append(lookup(pngsTable[id1][neededPosInds[pos]], pos, lookupSource))
                vec2.append(lookup(pngsTable[id2][neededPosInds[pos]], pos, lookupSource))
            total += weightedEuclidean(vec1, vec2, weights) / genDistMatrix[id1, id2]
        volatilities.append(total / len(genDistMatrix))
    else:
        for pos in neededPosInds:
            total = 0
            for id1, id2 in genDistMatrix:
                v1 = lookup(pngsTable[id1][neededPosInds[pos]], pos, lookupSource)
                v2 = lookup(pngsTable[id2][neededPosInds[pos]], pos, lookupSource)
                total += (v1 - v2) ** 2 / genDistMatrix[id1, id2]
            volatilities.append(total / len(genDistMatrix))
    return volatilities


# set up weight scale
if site == '2G12':
    weights = np.array([0.626, 0.81, 0.267, 0.45, 0.699])
elif site == '2F5':
    weights = np.array([1.242249, 0.86661, 0.219614, 0.855693, 0.789327])
elif site == 'all':
    pass
else:
    raise Exception("Unspecified calculation type")

# calculate all volatilities
pngsTable, allPosDict = readPngs(pngsTableFileName)  # read pngs table into memroy
neededPosInds = getNeededPositions(site, allPosDict)
for patientFileName in distanceFileNames:

    # read genetic distance matrix into memory
    genDistMatrix = readGeneticDistance(geneticDistanceFolderName + patientFileName)
    if len(genDistMatrix) == 0:  # empty hyphy output file due to too few envelopes in time point
        continue

    volatility = calcVolatility(genDistMatrix, pngsTable, neededPosInds, isFiveD, lookupSource)  # calculate volatility

    uniqueIds = set()  # calculate number of envelopes as number of unique identifiers in the distance matrix
    for id1, id2 in genDistMatrix:
        uniqueIds.add(id1)
        uniqueIds.add(id2)
    numEnvs = len(uniqueIds)
    if skipLessThanTwoVolatilities and numEnvs < 2:
        continue

    # generate output string id for the patient
    for key, _ in genDistMatrix:
        patientIdList = key.split('.')[:headerOffset - 1]
        break

    # write output to display
    strVolatility = [str(vol) for vol in volatility]
    print(" ".join([str(numEnvs)] + patientIdList + strVolatility))
