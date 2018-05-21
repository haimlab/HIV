import csv
import os
import re
import argparse
import itertools
import math
import sys

# filter positive integers used for bounds
appearedBounds = []
def bounds(string):
    if not string.isdigit():
        raise argparse.ArgumentTypeError('%s is not an integer' % string)
    val = int(string) 
    if val < 1:
        raise argparse.ArgumentTypeError('%s is not positive' % string)
    elif val in appearedBounds:
        raise argparse.ArgumentTypeError('%s is duplicated bound' % string)
    else: # check spacing between bounds
        for bound in appearedBounds:
            if abs(val - bound) < 2:
                raise argparse.ArgumentTypeError('0-sized bin')
    return val

# parse command line inputs
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bounds', nargs = '+', required = True, type = bounds)
parser.add_argument('-i', '--inFile', nargs = 1, required = True)
parser.add_argument('-o', '--outFile', nargs = '?', default = 'output.tab')
parser.add_argument('-f', '--falseVal', nargs = '?', default = 0.001, type = float)
try:
    args = parser.parse_args()
except:
    sys.exit('please double check input')

# helper method to process bin ranges
def parseBinRange(bounds):
    bounds.sort()
    curStart = 0
    binRanges = []
    for bound in bounds:
        binRanges.append((curStart, bound))
        curStart = bound + 1
    return binRanges

# get bin ranges
binRanges = parseBinRange(args.bounds)

# constants
dt = 1
posStart = 2
digitsPattern = re.compile('\d+')

# calculate number of combination pairs
def nCr(n, r):
    if type(n) != int or type(r) != int:
        raise Exception('cannot compute combination on non-integer')
    elif r > n:
        raise Exception('r must be no larger than n in combination')
    return math.factorial(n) / (math.factorial(r) * math.factorial(n - r))

# stores values of four quadrants used to calculate fisher exact test for a
# particular bin
class Position:
    def __init__(self, position):
        self.position = int(position)
        self.q1 = 0 # values of the four quadrants
        self.q2 = 0
        self.q3 = 0
        self.q4 = 0        
    def add(self, v1, v2): # add a pair of values (A, B) to this bin
        falseVal = args.falseVal
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

# represents a bin and all positions. The positions are created by considering
# only values that fall into this bin.
class Bin:
    def __init__(self, interval):
        self.start = interval[0] # interval of the bin
        self.end = interval[1]
        self.positions = {}
    def intervalContains(self, t):
        return t <= self.end and t >= self.start
    def toList(self): # turn the bin into a list representation, for csv writer
        myList = []
        keyList = []
        myList.append(str(self.start) + ' - ' + str(self.end))
        for pos in allPos:
            myList.append(self.positions[pos].calcFisher())
        return myList
    def add(self, position, val1, val2): # add a new value pair to this bin
        try:
            self.positions[position].add(val1, val2)
        except KeyError:
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
allPos = []
with open(args.inFile[0]) as inFile:
    reader = csv.reader(inFile)
    isFirstRow = True
    for row in reader:
        if isFirstRow:
            for i in range(0, len(row)):
                try:
                    pos = int(digitsPattern.match(row[i]).group(0))
                except AttributeError:
                    continue
                try:
                    posDict[pos].append(i)
                except KeyError:
                    posDict[pos] = [i]
                if not pos in allPos:
                    allPos.append(pos)
            allPos.sort()
            isFirstRow = False
            continue
        for b in bins:
            if b.intervalContains(int(row[dt])):
                for i in posDict:
                    b.add(i, float(row[posDict[i][0]]), float(row[posDict[i][1]]))
                break

# write output
with open(args.outFile, 'w') as outFile:
    writer = csv.writer(outFile, lineterminator = '\n', delimiter = '\t')
    header = ['bins']
    for pos in allPos:
        header.append(str(pos))
    writer.writerow(header)
    for b in bins:
        writer.writerow(b.toList())