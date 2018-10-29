chimport csv

# inputs
ranges = [-0.65, -0.5, -0.35, -0.25, 0.25, 0.35, 0.5, 0.65]  # range of inputs
# ranges = [-5.67067, -4.54392, -2.95569, -1.50762]
baseColor = '#0000000079e7'
colors = ['#30c330c3ffff', '#6db66db6ffff', '#9e799e79ffff', 'white', '#eebd861f861f', '#eebd43104310', '#eebd00000000', '#592b00000000']  # blue to red
# colors = ['#9e799e79ffff', '#55555555ffff', '#30c330c3ffff', '#0000000079e7']  # light to dark blue
# colors = ['#eebd861f861f', '#eebd43104310', '#eebd00000000', '#592b00000000']  # light to dark red
classificationValuesFileName = "..\\Colorings\\7.18.18_Delta_Volatility_And_Entropy\\values_for_classification.csv"
POS = 0  # column index of position in input file
DATA = 2  # column index of data to use for classification in input file
# end of inputs

# constants
trimerFileName = '~\\Desktop\\Trimer Coloring\\Colorings\\trimer_model.py'
commandFileName = 'tempCommand.cmd'
bashScriptFileName = 'tempBash.bat'
# end of constants

bins = [[] for i in range(0, len(ranges))]


# helper functions
def decideBin(position, halfLife):
    for i in range(0, len(ranges) - 1):
        if ranges[i] <= halfLife < ranges[i + 1]:
            bins[i].append(position)
            return True
    if halfLife > ranges[-1]:
        bins[-1].append(position)
        return True


# read through half life file and classify into bins
with open(classificationValuesFileName) as valuesFile:
    reader = csv.reader(valuesFile)
    firstLine = True
    for row in reader:
        if firstLine:
            firstLine = False
            continue
        if not row[DATA].strip() == "":
            decideBin(int(row[POS]), float(row[DATA]))

# classify positions that share the same color but with different numbers
POS1 = 5
POS2 = 13
moleculeFileNames = ["Other_Script_Constants_(Do_Not_Change)\\molecule" + str(i) + ".csv" for i in range(1, 4)]
for fileName in moleculeFileNames:
    with open(fileName) as file:
        reader = csv.reader(file)
        firstRow = True
        for row in reader:
            if firstRow:
                firstRow = False
                continue
            for b in bins:
                try:
                    p = int(row[POS1])
                except ValueError:
                    continue
                if p in b:
                    try:
                        b.append(int(row[POS2][1:]))
                    except ValueError:
                        continue

# decdie the ranges to write to out files
rangesToWrite = []
for i in range(0, len(ranges) - 1):
    rangesToWrite.append("(" + str(ranges[i]) + " - " + str(ranges[i + 1]) + ")")
rangesToWrite.append("(" + str(ranges[-1]) + "+)")

# write classified positions to a text file
# with open(outFileName, 'w') as outFile:
#    for i in range(0, len(bins)):
#         curBin = list(set(bins[i]))
#         outFile.write(rangesToWrite[i])
#         outFile.write(" ")
#        for j in range(0, len(curBin) - 1):
#            outFile.write(str(curBin[j]) + ", ")
#        outFile.write(str(curBin[-1]) + "\n")

# writing the according script
with open(bashScriptFileName, 'w') as bashScriptFile:
    bashScriptFile.write("set PATH=%PATH%;\"C:\\Program Files\\Chimera 1.12\\bin\"\n")
    bashScriptFile.write('chimera ' + commandFileName)


# write the according chimera command file
with open(commandFileName, 'w') as commandFile:
    commandFile.write('open ' + trimerFileName + '\n')
    commandFile.write('color ' + baseColor + '\n')
    for i in range(0, len(bins)):
        b = list(set(bins[i])) # remove duplicates and sort
        b.sort()
        bStr = [str(p) for p in b]
        commandFile.write('color ' + colors[i] + ' :' + ','.join(bStr) + '\n')

# execute the scripts we just created using python
from subprocess import Popen
p = Popen(bashScriptFileName)
stdout, stderr = p.communicate()

# remove temp files
import os
os.remove(commandFileName)
os.remove(bashScriptFileName)