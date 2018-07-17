import os
from look_up_tables import lookup
import math
import numpy as np

# all input variables
fiveD = True
vaccineSite = '2F5'  # or 2F5
geneticDistanceFolderName = "C:\\Users\\Rentian Dong\\Desktop\\AE_5D_And_Regular_Volatility\\Genetic_Distance" \
                            "\\Hyphy_Output\\"
pngsTableFileName = "C:\\Users\\Rentian Dong\\Desktop\\AE_5D_And_Regular_Volatility\\Z_Conversion\\converted.txt"
# end of all input variables

# get all distance file absolute path
distanceFileNames = os.listdir(geneticDistanceFolderName)


# takes in two vectors, and multiply number for each dimension
def weightedEuclidean(diffVector, weightVector):
    if len(diffVector) != len(weightVector):
        raise Exception("Dimensions of input vectors mismatch")
    mySum = 0
    for i in range(0, len(diffVector)):
        mySum += diffVector[i] ** 2 * weightVector[i]
    return mySum ** .5


# read table into memory
table = []
with open(pngsTableFileName, 'r') as tableFile:
    for line in tableFile:
        table.append(line)

# set up weight scale
if fiveD:
    if vaccineSite == '2G12':
        weights = np.array([0.626, 0.81, 0.267, 0.45, 0.699])
    elif vaccineSite == '2F5':
        weights = np.array([1.242249, 0.86661, 0.219614, 0.855693, 0.789327])
    else:
        raise Exception("Unspecified calculation type")


# mysterious function that does something special if a patient name starts with 2e
def mystic2eHandler(patientName):
    patientName = patientName[2:]  # remove 2e from patient name
    flag = False
    processedString = ''
    for char in patientName:
        if char == '-':
            flag = True
        elif flag:
            if char == '.':
                flag = False
                processedString += '.'
        else:
            processedString += char
    return processedString


# HyPhy Philips Style Matrix has number of Sequences listed in the first position of output
# This code brings the matrix into python
# read genetic distance table as a matrix into memory
# returns (allIdentifiders, distanceMatrix, #envelopes)
def readGeneticDistance(distanceFileName):
    identifiers = set()
    geneticDistanceMatrix = {}

    with open(distanceFileName, 'r') as distanceFile:
        isFirstRow = True
        for row in distanceFile:
            if isFirstRow:  # skip first row
                isFirstRow = False
                continue

            t = row.split()

            if nam[:2] == '2e':
                z0 = mystic2eHandler(t[0])
                z1 = mystic2eHandler(t[1])
            else:
                z0 = t[0].replace('-', '')
                z1 = t[1].replace('-', '')

            if z0 not in identifiers:
                identifiers.add(z0)
            if z1 not in identifiers:
                identifiers.add(z1)

            genDist = float(t[2])
            if genDist <= 0:
                raise Exception('non positive genetic distance')
            geneticDistanceMatrix[(z0, z1)] = genDist
            geneticDistanceMatrix[(z1, z0)] = genDist

        numEnvelopes = len(identifiers)
        return identifiers, geneticDistanceMatrix, numEnvelopes


# read in the pngs file rows whose identifier mathces that provided
def readPngs(allIds):
    # read in pngs table
    lines = []  # list of linse of pngs table, where each line separated in to list as space delimited
    years = []
    phases = []

    isFirstRow = True
    for line in table:
        if len(line) == 0:
            continue
        t = line.split()
        if isFirstRow:
            isFirstRow = False
            if t[0] != 'Country':
                raise Exception('header format seems incorrect')
            lines.append(t)
            continue

        ID = t[0] + '.' + t[1] + '.' + t[2].split("-")[0] + '.' + t[3]
        if ID not in allIds:
            continue
        if t[1] not in years:
            years.append(t[1])
        try:
            if t[2].split("-")[1] not in phases:
                phases.append(t[2].split("-")[1])
        except IndexError:
            if "" not in phases:
                phases.append("")
        lines.append(t)

    return lines, years, phases


# iterate and calculate all volatilities
for nam in distanceFileNames:

    # read distance matrix and pngs table
    allIds, genDistMatrix, numEnvelopes = readGeneticDistance(geneticDistanceFolderName + nam)
    if numEnvelopes < 2:  # skip samples with 2 or less envelopes
        continue
    lines, years, phases = readPngs(allIds)
    positions = lines[0][4:]  # all positions present in the pngs file

    if vaccineSite == '2G12':  # positions needed for this instance of calculation
        neededPositions = {295, 332, 339, 392, 448}
        posIndices = dict()
        for i in range(0, len(positions)):
            if int(positions[i]) in neededPositions:
                posIndices[int(positions[i])] = i + 4
    elif vaccineSite == '2F5':
        neededPositions = {662, 663, 664, 665, 667}
        posIndices = {}
        for i in range(0, len(positions)):
            if int(positions[i]) in neededPositions:
                posIndices[int(positions[i])] = i + 4
    elif vaccineSite == 'Hydropathy':  # every position present in the file
        neededPositions = set([int(pos) for pos in positions])
    else:
        raise Exception("unidentified position set")

    num = 0
    vols = []
    for year in years:
        for phase in phases:
            for pos in positions:
                pos = int(pos)
                if pos not in neededPositions:  # calculate only the specified vaccine sites if required
                    continue
                envIds = []
                PatNum = [0] * numEnvelopes
                positionVals = []
                count = 0
                yr_counter = 0
                for line in lines:
                    skip = False
                    if len(line) == 0:
                        break
                    try:
                        if line[1] != year or line[2].split("-")[1] != phase:
                            yr_counter += 1
                            skip = True
                    except IndexError:
                        pass

                    if count < 1:  # ignoring header row of pngs table
                        count += 1
                        continue
                    envIds.append(line[0] + '.' + line[1] + '.' + line[2].split("-")[0] + '.' + line[3])
                    if skip:
                        positionVals.append(None)
                    elif fiveD:
                        arr = np.array([lookup(line[posIndices[key]], key, vaccineSite) for key in posIndices])
                        positionVals.append(arr)
                    else:
                        positionVals.append(lookup(line[posIndices[pos]], pos, vaccineSite))
                    count += 1

                totalVol = 0
                numCalcs = 0
                for i in range(0, numEnvelopes):
                    for j in range(0, i):
                        if positionVals[i] is None or positionVals[j] is None:
                            continue
                        diffs = positionVals[i] - positionVals[j]
                        if fiveD:
                            phenotypicDist = weightedEuclidean(diffs, weights)
                        else:
                            phenotypicDist = math.pow(diffs, 2)
                        totalVol += phenotypicDist / genDistMatrix[(envIds[i], envIds[j])]
                        numCalcs += 1
                vols.append(str(totalVol / numCalcs))

                num = count - yr_counter

                if fiveD:
                    break

            vol_string = " ".join(vols)
            try:
                if num > 1:
                    resultStrList = [str(num), line[0], year, line[2].split("-")[0], phase, line[3]]
            except KeyError:
                resultStrList = [str(num), '-', year, nam, phase, '-', vols[0]]
            if fiveD:
                resultStrList.append(vols[0])
            else:
                resultStrList.append(vol_string)
            print(' '.join(resultStrList))
