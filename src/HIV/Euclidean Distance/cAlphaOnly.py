import csv
import re

# (CONSTANTS) indices of elements
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

# helper method for computing euclidean distance
def euclidean3D(v1,v2):
    if len(v1) != 3 or len(v2) != 3:
        raise Exception("vector lengths incorrect")
    mySum = 0
    for i in range(0, 3):
        mySum += (v1[i] - v2[i]) ** 2
    return mySum ** 0.5

def genKey(p1, p2):
    if int(p1) > int(p2):
        return p2 + "-" + p1
    else:
        return p1 + "-" + p2

files = []
for p in filePaths:
    file = open(p, "r")
    files.append(file)

readers = []
for f in files:
    reader = csv.reader(f)
    readers.append(reader)

# create dictionaries for c-alphas of each molecule
cAlphas = {}
header = []
for r in readers:
    # extract all c-alpha coordinates
    for row in r:
        if row[COL1] == "ATOM" and row[ATOMTYPE] == "CA":
            chain = row[CHAIN]
            pos = row[POS]
            if pos.isdigit():
                if not int(pos) in header:
                    header.append(int(pos))
            key = chain + "," + pos
            cAlphas[key] = [float(row[X]), float(row[Y]), float(row[Z])]
header.sort()
temp = []
for element in header:
    temp.append(str(element))
header = temp

# calculate distance of all c-alphas
results = {}
for base in cAlphas:
    for other in cAlphas:
        basePair = base.split(",")
        otherPair = other.split(",")
        if len(basePair) != 2 or len(otherPair) != 2:
            raise Excpetion("incorrect key format!")
        pos1 = basePair[1]
        pos2 = otherPair[1]
        # skip comparing to oneself
        if pos1 == pos2 and basePair[0] == otherPair[0]:
            continue
        # skip insertions
        if (not pos1.isdigit()) or (not pos2.isdigit()):
            continue
        int1 = int(pos1)
        int2 = int(pos2)
        if int1 < int2:
            resultKey = pos1 + "-" + pos2
        else:
            resultKey = pos2 + "-" + pos1
        d = euclidean3D(cAlphas[base], cAlphas[other])
        try:
            oldVal = results[resultKey]
            if d < oldVal:
                results[resultKey] = d
        except KeyError:
            results[resultKey] = d

# write results to output
with open("out.csv", "w") as out:
    out.write(",")
    for pos in header:
        out.write(pos + ",")
    out.write("\n")
    # rows
    for rowPos in header:
        out.write(rowPos + ",")
        # pos on a row
        for pos in header:
            if int(pos) <= int(rowPos):
                out.write(",")
            else:
                key = genKey(rowPos, pos)
                out.write(str(results[key]))
                out.write(",")
        out.write("\n")
        
for f in files:
    f.close()