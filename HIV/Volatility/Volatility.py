import os
from look_up_tables import lookup
import math
import numpy as np

# all input variables
fiveD = True
vaccineSite = '2F5' # or 2F5
geneticDistanceFolderName = "C:\\Users\\rdong6\\Desktop\\AE_5D_And_Regular_Volatility\\Genetic_Distance\\Hyphy_Output\\"
pngsTableFileName = "C:\\Users\\rdong6\\Desktop\\AE_5D_And_Regular_Volatility\\Z_Conversion\\converted.txt"
# end of all input variables

# get all distance file absolute path
distanceFileNames = os.listdir(geneticDistanceFolderName)

# takes in two vectors, and multiply number for each dimension
def weightedEuclidean(diffVector, weights):
    if len(diffVector) != len(weights):
        raise Exception("Dimensions of input vectors mismatch")
    sum = 0
    for i in range(0, len(diffVector)):
        sum += diffVector[i] ** 2 * weights[i]
    return sum ** .5

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

# mysterious function that does something specifial if a patient name starts with 2e
def mystic2eHandler(patientName):
    patientName = patientName[2:] # remove 2e from patient name
    flag = False
    processedString = ''
    for char in t[0]:
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
            if isFirstRow: # skip first row
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

# iterate and calculate all volatilities
for nam in distanceFileNames:

    # read in distance matrix and its identifiers
    allIds, genDistMatrix, numEnvelopes = readGeneticDistance(geneticDistanceFolderName + nam)
    if numEnvelopes < 2: # skip samples with 2 or less envelopes
        continue

    # read in pngs table
    lines = [] # list of linse of pngs table, where each line separated in to list as space delimited
    years = []
    phases = []

    isFirstRow = True
    for line in table:
        if len(line) == 0:
            break
        t = line.split()
        if isFirstRow:
            isFirstRow = False
            if t[0] != 'Country':
                raise Exception('header format seems incorrect')
            lines.append(t)
            continue

        id = t[0] + '.' + t[1] + '.' + t[2].split("-")[0] + '.' + t[3]
        if id not in allIds:
            continue
        if t[1] not in years:
            years.append(t[1])
        try:
            if t[2].split("-")[1] not in phases:
                phases.append(t[2].split("-")[1])
        except:
            if "" not in phases:
                phases.append("")
        lines.append(t)

    positions = lines[0][4:]  # all positions present in the pngs file
    if vaccineSite == '2G12': # positions needed for this instance of calculation
        neededPositions = {295, 332, 339, 392, 448}
        neededPositionsIndex = dict()
        for i in range(0, len(positions)):
            if int(positions[i]) in neededPositions:
                neededPositionsIndex[int(positions[i])] = i + 4
    elif vaccineSite == '2F5':
        neededPositions = {662, 663, 664, 665, 667}
        neededPositionsIndex = {}
        for i in range(0, len(positions)):
            if int(positions[i]) in neededPositions:
                neededPositionsIndex[int(positions[i])] = i + 4
    elif vaccineSite == 'Hydropathy': # every position present in the file
        neededPositions = set([int(pos) for pos in positions])
    else:
        raise Exception("unidentified position set")

    num = 0
    vols = [0] * len(neededPositions)

    for year in years:
        for phase in phases:
            posCount = 0
            for pos in positions:
                # calculate only the specified vaccine sites if required
                if int(pos) not in neededPositions:
                    continue
                envIds = [0] * numEnvelopes
                PatNum = [0] * numEnvelopes
                positionVals = [0] * numEnvelopes
                count = 0
                yr_counter = 0
                for line in lines:
                    skipper = False
                    already_skip = True
                    if len(line) == 0:
                        break
                    if line[1] != year:
                        yr_counter += 1
                        skipper = True
                        already_skip = False
                    try:
                        if line[2].split("-")[1] != phase:
                            if already_skip:
                                yr_counter += 1
                                skipper = True
                    except IndexError:
                        pass

                    if count < 1:
                        count += 1
                        continue
                    envIds[count - 1] = line[0] + '.' + line[1] + '.' + line[2].split("-")[0] + '.' + line[3]
                    if(len(line) > (posCount+4)):
                        positionVals[count - 1] = lookup(line[posCount + 4], int(pos), vaccineSite) # hydropathy
                    if fiveD:
                        #print(envIds[count - 1], flush=True)
                        #string = ''
                        #for key in neededPositionsIndex:
                        #    string += line[neededPositionsIndex[key]]
                        #print(string, flush=True)
                        positionVals[count - 1] = np.array([
                            lookup(line[neededPositionsIndex[key]], key, vaccineSite) for key in neededPositionsIndex
                        ])
                        #print(positionVals[count - 1], flush=True)
                        #print('', flush=True)

                    if skipper:
                        if fiveD:
                            positionVals[count - 1] = np.array([-1000])
                        else:
                            positionVals[count - 1] = -1000
                    count += 1

                deltaVals = [[0.0 for _ in range(numEnvelopes)] for _ in range(numEnvelopes)]
                volVals = [[0.0 for _ in range(numEnvelopes)] for _ in range(numEnvelopes)]
                totalVol = 0
                numCalcs = 0
                for i in range(0, numEnvelopes):
                    for j in range(0, numEnvelopes):
                        if i == j:
                            break
                        if fiveD:
                            try:
                                if positionVals[i][0] == -1000 or positionVals[j][0] == -1000:
                                    continue
                            except:
                                pass
                        else:
                            if positionVals[i] == -1000 or positionVals[j] == -1000:
                                continue

                        diffs = positionVals[i] - positionVals[j]
                        if fiveD:
                            phenotypicDist = weightedEuclidean(diffs, weights)
                        else:
                            phenotypicDist = math.pow(diffs, 2)
                        totalVol += phenotypicDist / genDistMatrix[(envIds[i], envIds[j])]
                        numCalcs += 1

                vols[posCount] = str(totalVol / numCalcs)
                posCount += 1

                num = count - yr_counter

                if fiveD:
                    break
            vol_string = " "
            for vol in vols:
                vol_string += str(vol) + ' '

            try:
                if num > 1:
                    if fiveD:
                        print(str(num) + ' ' + line[0] + ' ' + year + ' ' + line[2].split("-")[0] + ' ' + phase + ' ' + line[3] + '  ' + vols[0])
                    else:
                        print(str(num) + ' ' + line[0] + ' ' + year + ' ' + line[2].split("-")[0] + ' ' + phase + ' ' + line[3] + ' ' + vol_string)
            except:
                if fiveD:
                    print(str(num) + " " + '-' + ' ' + year + ' ' + nam + ' ' + phase + ' ' + '-' + '  ' + vols[0])
                else:
                    print(str(num) +" "+ '-' + ' ' + year +' '+ nam +' '+ phase+' '+'-'+' '+vol_string)