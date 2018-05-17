import csv
import re
import sys

# indicies of elements
COL1 = 0
ATOMTYPE = 2
AA = 3
CHAIN = 4
POS = 5
X = 6
Y = 7
Z = 8

# path to molecules 1, 2 and 3
filePaths = ["m1.csv", # BG
                 "m2.csv", # CH
                 "m3.csv"] # DI

# DS to store coordinates of all atoms at one position
class Position:
    def __init__(self, p):
        self.position = p
        self.coords = []
    def addCoord(self, toAdd):
        if not len(toAdd) == 3:
            raise Exception("invalid coordinates")
        self.coords.append(toAdd)

# take two position instance as input, and output closest distance as float
def findClosest(p1, p2):
    minimum = sys.float_info.max
    for c1 in p1.coords:
        for c2 in p2.coords:
            newDist = euclidean3D(c1, c2)
            if (newDist < minimum):
                minimum = newDist
    return minimum

# helper method for computing euclidean distance
def euclidean3D(v1,v2):
    if len(v1) != 3 or len(v2) != 3:
        raise Exception("vector lengths incorrect")
    mySum = 0
    for i in range(0, 3):
        mySum += (v1[i] - v2[i]) ** 2
    return mySum ** 0.5

def genKey(p1, p2):
    if p1 > p2:
        return str(p2) + "-" + str(p1)
    else:
        return str(p1) + "-" + str(p2)

files = []
for p in filePaths:
    file = open(p, "r")
    files.append(file)

readers = []
for f in files:
    reader = csv.reader(f)
    readers.append(reader)

# create dictionaries for c-alphas of each molecule
allPosCoords = {}
header = []
for r in readers:
    # extract all c-alpha coordinates
    for row in r:
        if row[COL1] == "ATOM" and row[AA] != "HOH" and row[AA] != "SO4":
            pos = row[POS]
            if pos.isdigit():
                pos = int(pos)
                if not pos in header:
                    header.append(pos)
                # if no DS for the position
                if not pos in allPosCoords:
                    newPos = Position(pos)
                    newPos.addCoord([float(row[X]), float(row[Y]), float(row[Z])])
                    allPosCoords[pos] = newPos
                # if DS for position already created
                else:
                    allPosCoords[pos].addCoord([float(row[X]), float(row[Y]), float(row[Z])])
header.sort()
print("finished parsing data")

# calculate distance of all c-alphas
results = {}
counter = 0
for pos1 in allPosCoords:
    for pos2 in allPosCoords:
        resultKey = genKey(pos1, pos2)
        print("calculating " + resultKey, flush = True)
        counter += 1
        print(counter)
        d = findClosest(allPosCoords[pos1],  allPosCoords[pos2])
        # print("calculated pair " + resultKey + " to be " + str(d), flush = True)
        try:
            oldVal = results[resultKey]
            if d < oldVal:
                results[resultKey] = d
        except KeyError:
            results[resultKey] = d
print("finished calculating distances")

# write results to output
with open("out.csv", "w") as out:
    out.write(",")
    for pos in header:
        out.write(str(pos) + ",")
    out.write("\n")
    # rows
    print(header)
    for rowPos in header:
        print("writing row " + str(rowPos), flush = True)
        # debug
        out.write(str(rowPos) + ",")
        # pos on a row
        for pos in header:
            if pos <= rowPos:
                out.write(",")
            else:
                key = genKey(rowPos, pos)
                out.write(str(results[key]))
                out.write(",")
        out.write("\n")
print("finished writing")
        
for f in files:
    f.close()