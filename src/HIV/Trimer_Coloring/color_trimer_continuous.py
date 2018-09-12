import csv
import colorsys
from subprocess import Popen
import os

# inputs
# [(minVal, maxVal, (h, s, v,), (h, s, v,))]

#ranges = [ # delta volatility
#   (-.65, -.25, (.681, 1, .6), (.681, 0.2, 1)),
#   (-.25, .25, (1, 0, 1), (1, 0, 1)), # always white
#   (.25, .65, (1, 0.2, 1), (1, 1, .6))
#]

#ranges = [ # delta entropy
#   (-.65, -.25, (.681, 1, .6), (.681, 0, 1)),
#   (-.25, .25, (1, 0, 1), (1, 0, 1)), # always white
#   (.25, .65, (1, 0, 1), (1, 1, .6))
#]

#ranges = [ # p6 clade B entropy
#    (0.217, 2.266, (.681, 0, 1), (.681, 1, .6))
#]

#ranges = [ # p6 clade B volatility
#    (-5.70853, -2.08034, (.681, 0, 1), (.681, 1, .6))
#]

#ranges = [ # p6 clade C entorpy
#    (0.218, 2.279, (1, 0, 1), (1, 1, .6))
#]

ranges = [ # p6 clade C volatility
    (-5.67041, -1.506, (1, 0, 1), (1, 1, .6))
]

classificationValuesFileName = "..\\values_for_classification.csv"
POS = 0  # column index of position in input file
DATA = 6  # column index of data to use for classification in input file
# end of inputs


# constants
POS_1600_FILE_NAMES = ["Other_Script_Constants_(Do_Not_Change)\\molecule" + str(i) + ".csv" for i in range(1, 4)]
TRIMER_FILE_NAME = 'Other_Script_Constants_(Do_Not_Change)\\trimer_model.py'

commandFileName = 'tempCommand.cmd'
batchScriptFileName = 'tempBatch.bat'
POS_INDEX = 5
POS1600_INDEX = 13
# end of constants

class ColorRange:
    def __init__(self, minVal, maxVal, minColor, maxColor):
        self.minVal = minVal
        self.maxVal = maxVal
        self.minColor = minColor
        self.maxColor = maxColor


# parse input color ranges
def parseColorRanges(ranges):
    colorRangeObjects = []
    for r in ranges:
        colorRangeObjects.append(ColorRange(r[0], r[1], r[2], r[3]))
    return colorRangeObjects


# linear interpolation
def linIntp(x1, x2, y1, y2, x):
    return (x - x1) * (y2 - y1) / (x2 - x1) + y1


# helper functions
def interpolateColor(minVal, maxVal, minColor, maxColor, val):
    newColor = []
    for i in range(0, 3):  # h, s, and v
        newColor.append(linIntp(minVal, maxVal, minColor[i], maxColor[i], val))
    return tuple(newColor)


def decideColor(value, colorRanges):
    # color belongs to one of the ranges
    for i in range(0, len(colorRanges)):
        if colorRanges[i].minVal <= value < colorRanges[i].maxVal:
            return interpolateColor(colorRanges[i].minVal, colorRanges[i].maxVal, colorRanges[i].minColor,
                                    colorRanges[i].maxColor, value)
    if value >= colorRanges[-1].maxVal:
        return colorRanges[-1].maxColor
    elif value < colorRanges[0].minVal:
        return colorRanges[0].minColor
    else:
        raise Exception("Unable to determine color for value: " + str(value))


# read input csv file into position -> value dictionary
def parsePosValFile(fileName, posCol, valCol):
    posValDict = {}
    with open(fileName, 'r') as file:
        reader = csv.reader(file)
        firstRow = True
        for row in reader:
            if firstRow:
                firstRow = False
                continue
            posValDict[int(row[posCol])] = float(row[valCol])
    return posValDict


# take in pos -> val dictionary and classify into  position -> RGB
def determineAllColors(posValDict, posTo1600Dict, colorRanges):
    posColorDict = {}
    for pos in posValDict:
        posColorDict[pos] = decideColor(posValDict[pos], colorRanges)
    posColorDict = addd1600Pos(posTo1600Dict, posColorDict)
    return posColorDict


# return dictionary of pos -> {pos} (set) mapping 1600s positions that should share the same color
def parse1600PosFiles(fileNames):
    posTo1600Dict = {}
    for fileName in fileNames:
        with open(fileName) as file:
            reader = csv.reader(file)
            firstRow = True
            for row in reader:
                if firstRow:
                    firstRow = False
                    continue
                try:  # extract position pair
                    pos = int(row[POS_INDEX])
                    pos1600 = int(row[POS1600_INDEX][1:])
                except ValueError:
                    continue
                try:  # put position pair into dictionary
                    posTo1600Dict[pos].add(pos1600)
                except KeyError:
                    posTo1600Dict[pos] = {pos1600}
    return posTo1600Dict


# add 1600 positions to pos -> color dictionary
def addd1600Pos(dict1600, dictColor):
    for pos in dict1600:
        for pos1600 in dict1600[pos]:
            dictColor[pos1600] = dictColor[pos]
    return dictColor


# write the batch script file
def writeBatchScript():
    with open(batchScriptFileName, 'w') as bashScriptFile:
        bashScriptFile.write("set PATH=%PATH%;\"C:\\Program Files\\Chimera 1.13\\bin\"\n")
        bashScriptFile.write('chimera ' + commandFileName)


# turn (h, s, v) into r, g, b as a string
def formatColorCode(hsvCode):
    rgbCode = colorsys.hsv_to_rgb(hsvCode[0], hsvCode[1], hsvCode[2])
    rgbStrCode = [str(val) for val in rgbCode]
    return ','.join(rgbStrCode)


# write the command file based on position color pairs
def writeCommandFile(posColorDict):
    with open(commandFileName, 'w') as commandFile:
        commandFile.write('open ' + TRIMER_FILE_NAME + '\n')
        commandFile.write('color green\n')
        for pos in posColorDict:
            colorStr = formatColorCode(posColorDict[pos])
            commandFile.write('color ' + colorStr + ' :' + str(pos) + '\n')


if __name__ == '__main__':

    # calculate the colors
    colorRanges = parseColorRanges(ranges)
    posToVal = parsePosValFile(classificationValuesFileName, POS, DATA)
    posTo1600Dict = parse1600PosFiles(POS_1600_FILE_NAMES)
    posToColor = determineAllColors(posToVal, posTo1600Dict, colorRanges)

    # write calculated colors to batch scripts excutable by chimera
    writeBatchScript()
    writeCommandFile(posToColor)

    # execute the scripts we just created
    p = Popen(batchScriptFileName)
    stdout, stderr = p.communicate()

    # delete the temporary scripts
    os.remove(commandFileName)
    os.remove(batchScriptFileName)
