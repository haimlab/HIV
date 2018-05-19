# fisher z-transformation
import csv
import os
import re

# input parameters
binRanges = [(1, 31), (32, 90), (91, 180), (101, 365), (366, 710), (730, 4845)] # both ends inclusive
folderPath = '.'
inFileName = 'testInput2.csv'
outFileName = 'out.csv'
falseVal = 0.001
posNeeded = ["295","332","339","392","386","448","412","413","137", \
                   "301","327","363","662","663","664","665","666","667", \
                   "668","669","670","671","672","673","674","675","676", \
                   "677","678","679","680","681","682","683","156","160", \
                   "165","167","168","169","170","171","173"]

# constants
dt = 1
posStart = 2
digitsPattern = re.compile('\d+')

def factorial(n):
    if type(n) != int:
        raise Exception('cannot compute factorial on non-integer')
    if n < 0:
        raise Exception('cannot compute factorial for negative number')
    if n == 1 or n == 0:
        return 1
    else:
        return n * factorial(n - 1)

def nCr(n, r):
    if type(n) != int or type(r) != int:
        raise Exception('cannot compute combination on non-integer')
    elif r > n:
        raise Exception('r must be no larger than n in combination')
    return factorial(n) / (factorial(r) * factorial(n - r))

class Position:
    def __init__(self, position):
        self.position = int(position)
        self.q1 = 0 # values of the four quadrants
        self.q2 = 0
        self.q3 = 0
        self.q4 = 0        
    def add(self, v1, v2): # add a pair of values (A, B) to this bin
        if v1 == falseVal and v2 == falseVal:
            self.q4 += 1
        elif v1 == falseVal and v2 != falseVal:
            self.q1 += 1
        elif v1 != falseVal and v2 == falseVal:
            self.q3 += 1
        else:
            self.q2 += 1
    def calcFisher(self): # calculate fisher using qudrant values of bin
        total = self.q1 + self.q2 + self.q3 + self.q4
        leftNumerator = nCr(self.q2 + self.q1, self.q2)
        rightNumerator = nCr(self.q3 + self.q4, self.q3)
        denominator = nCr(total, self.q2 + self.q3)
        fisherP = leftNumerator * rightNumerator / denominator
        return fisherP

class Bin:
    def __init__(self, interval):
        self.start = interval[0] # interval of the bin
        self.end = interval[1]
        self.positions = {}
    def intervalContains(self, t):
        return t <= self.end and t >= self.start
    def toList(self):
        myList = []
        myList.append(str(self.start) + ' - ' + str(self.end))
        for pos in posNeeded:
            pos = int(pos)
            myList.append(self.positions[pos].calcFisher())
        return myList
    def add(self, position, val1, val2):
        print('adding to positions' + str(position), flush = True)
        try:
            print('position already exist', flush = True)
            self.positions[position].add(val1, val2)
        except KeyError:
            print('creating a new position', flush = True)
            newPos = Position(position)
            newPos.add(val1, val2)
            self.positions[position] = newPos

# initialize variables to hold the data
bins = [] # this should and will remain sorted
for binRange in binRanges:
    binRange = (int(binRange[0]), int(binRange[1]))
    bins.append(Bin(binRange))

# parse input
posDict = {}
inFilePath = os.path.join(folderPath, inFileName)
inFile = open(inFilePath, 'r')
reader = csv.reader(inFile)
isFirstRow = True
for row in reader:
    if isFirstRow:
        for i in range(0, len(row)):
            try:
                pos = digitsPattern.match(row[i]).group(0)
            except AttributeError:
                continue
            if pos in posNeeded:
                try:
                    posDict[int(pos)].append(i)
                except KeyError:
                    posDict[int(pos)] = [i]
        isFirstRow = False
        continue
    print('checking row with dt: ' + row[dt])
    for b in bins:
        if b.intervalContains(int(row[dt])):
            print('contained in: ' + str(b.start) + ' - ' + str(b.end))
            for i in posDict:
                b.add(i, float(row[posDict[i][0]]), float(row[posDict[i][1]]))
            break

# write output
outFilePath = os.path.join(folderPath, outFileName)
outFile = open(outFilePath, 'w')
writer = csv.writer(outFile, lineterminator = '\n')
header = ['bins']
for pos in posNeeded:
    header.append(str(pos))
writer.writerow(header)
for b in bins:
    writer.writerow(b.toList())

# close files
inFile.close()
outFile.close()