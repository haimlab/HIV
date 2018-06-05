import random
import csv

# input parameters
inFileName = 'decideBins.txt'

# constants
POS_START = 2
PATIENT_INDEX = 0
DT_INDEX = 1
TIME_THRESH = 380

# to keep track of # of each transmission and the patient in each file
class FileStatus:
    def __init__(self):
        self.t1 = 0 # vt = 0, vt+1 = 0
        self.t2 = 0 # vt = 0, vt+1 = 1
        self.t3 = 0 # vt = 1, vt+1 = 0
        self.t4 = 0 # vt = 1, vt+1 = 1
        self.patientList = [] # list of patient names alloted to this file
    def combine(self, otherStatus):
        newStatus = FileStatus()
        newStatus.t1 = self.t1 + otherStatus.t1
        newStatus.t2 = self.t2 + otherStatus.t2
        newStatus.t3 = self.t3 + otherStatus.t3
        newStatus.t4 = self.t4 + otherStatus.t4
        newStatus.patientList = self.patientList + otherStatus.patientList
        return newStatus
    def getSum(self):
        return self.t1 + self.t2 + self.t3 + self.t4
    def __str__(self):
        return str(self.getSum())

def checkDuplicate(list1, list2):
    for element in list1:
        if element in list2:
            return True

# takes a header as a viewer and returns 0 based indicies of needed fields,
# generates an error if a particular field is missing
def getIndices(header, fields):
    indices = {}
    count = 0
    for element in header:
        if element in fields:
            try:
                if indices[element]:
                    raise Exception('duplicate field names')
            except KeyError:
                indices[count] = element
        count += 1
    return indices

# takes a header as a viewer and returns the index of first position
def findPosStart(header):
    count = 0
    for element in header:
        if element == 1:
            return count
        count += 1
    
# calculate the difference of the four quadrants    
def calcDiff(s1, s2):
    total = 0
    total += abs(s1.t1 - s2.t1)
    total += abs(s1.t2 - s2.t2)
    total += abs(s1.t3 - s2.t3)
    total += abs(s1.t4 - s2.t4)
    return total

def decideFile(s1, s2, p):
    newS1 = s1.combine(p)
    newS2 = s2.combine(p)
    return newS1.getSum() - newS2.getSum() < 0

# remove the leading and trailing '' elements in a list
def stripList(myList):
    while len(myList) > 0:
        if myList[0] == '':
            myList = myList[1:]
        else:
            break
    while len(myList) > 0:
        if myList[-1] == '':
            myList = myList[:len(myList) - 1]
        else:
            break
    return myList
    
# update a file status with the given row
def upDateStatus(status, row):
    # ignore if dT is larger than a threshold (currently 380)
    if int(row[DT_INDEX]) > 380:
        return
    row = stripList(row)
    for i in range(POS_START, len(row), 2):
        if row[i] == '0.001':
            v2 = 0
        else:
            v2 = float(row[i])
        if row[i + 1] == '0.001':
            v1 = 0
        else:
            v1 = float(row[i + 1])
        if v1 == 0 and v2 == 0:
            status.t1 += 1
        elif v1 == 0 and v2 > 0:
            status.t2 += 1
        elif v1 > 0 and v2 == 0:
            status.t3 += 1
        elif v1 > 0 and v2 > 0:
            status.t4 += 1
        else:
            raise Exception("Uunidentified Status")
    if row[PATIENT_INDEX] not in status.patientList:
        status.patientList.append(row[PATIENT_INDEX])

def writeStatus(status1, status2, preStatus):
    if decideFile(status1, status2, preStatus):
        status1 = status1.combine(preStatus)
    else:
        status2 = status2.combine(preStatus)
    prePatient = line[PATIENT_INDEX]
    preStatus = FileStatus()
    return [status1, status2, preStatus]

# assumes 'row' = first row, output the column index with patient names
def getPatientIndex(row):
    for i in range(0, len(row)):
        if line[i] == 'Patient':
            return i   

# initialize DS holding status of test and validation set
status1 = FileStatus()
status2 = FileStatus()
prePatient = ''
preStatus = FileStatus()

# read input file and determine which patient goes to which file
isFirstRow = True
with open(inFileName, 'r') as inFile:
    reader = csv.reader(inFile)
    for line in reader:
        if isFirstRow:
            # overwrite if found patient column
            newInd = getPatientIndex(line)
            if newInd != None:
                PATIENT_INDEX = newInd
            isFirstRow = False
            continue
        if line[PATIENT_INDEX] != prePatient:
            if not prePatient == '':
                [status1, status2, preStatus] = writeStatus(status1, status2, preStatus)
            prePatient = line[PATIENT_INDEX]
        # update current status
        upDateStatus(preStatus, line)
    # capture the last patient
    [status1, status2, preStatus] = writeStatus(status1, status2, preStatus)

# read input file and then split
with open('randFile1.csv', 'w') as file1, open('randFile2.csv', 'w') as file2, open(inFileName, 'r') as inFile:
    reader = csv.reader(inFile)
    writer1 = csv.writer(file1, lineterminator = '\n')
    writer2 = csv.writer(file2, lineterminator = '\n')
    isFirstRow = True
    for line in reader:
        # writer first row to both output files
        if isFirstRow:
            writer1.writerow(line)
            writer2.writerow(line)
            isFirstRow = False
            continue
        if int(line[DT_INDEX]) < 380:
            if line[PATIENT_INDEX] in status1.patientList:
                writer1.writerow(line)
            else:
                writer2.writerow(line)