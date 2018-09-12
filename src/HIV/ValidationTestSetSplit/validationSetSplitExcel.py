import pyexcel
import random
import csv

# input parameters
fields = ['Country', 'Year', 'Patient', 'Accession']
inFile = 'test.xlsx'

# to keep track of # of each transmission and the patient in each file
class FileStatus:
    def __init__(self):
        self.t1 = 0 # vt = 0, vt+1 = 0
        self.t2 = 0 # vt = 0, vt+1 = 1
        self.t3 = 0 # vt = 1, vt+1 = 0
        self.t4 = 0 # vt = 1, vt+1 = 1
        self.patientList = [] # list of patient names alloted to this file
        
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

# initialize DS holding status of test and validation set
status1 = FileStatus
status2 = FileStatus

# read input file and then split
book = pyexcel.get_book(file_name = 'test.xlsx')
with open('randFile1.csv', 'w') as file1, open('randFile2.csv', 'w') as file2:
    writer1 = csv.writer(file1)
    writer2 = csv.writer(file2)
    for sheet in book:
        isFirstRow = True
        for row in sheet:
            # handle header
            if isFirstRow:
                isFirstRow = False
                indices = getIndices(row, fields)
                posStart = findPosStart(row)
            
            else:
                # initialize fields
                newRow = []
                count = 0
                newHead = {}
                
                # get values of the required fields
                for cell in row:
                    if count in indices:
                        newHead[indices[count]] = cell
                    if len(newHead) == len(fields):
                        break;
                    count += 1
                
                # append values of required fields to new row, in fixed order
                for fieldName in fields:
                    newRow.append(newHead[fieldName])
                
                # put values under positions into new row
                count = 0
                for cell in row:
                    if count > posStart - 1:
                        newRow.append(cell)
                    count += 1
                
                # decide which file shoud this row go to randomly
                if newHead['Patient'] in status1.patientList: # ensure same patient goes to same file
                    writer1.writerrow(newRow)
                elif newHead['Patient'] in status2.patientList:
                    writer2.writerrow(newRow)
                else:
                    if random.randint(0, 1): # TODO, need to satisfy requirements of validation and test sets
                        writer1.writerow(newRow) # TODO, YOU ARE NOT UPDATING PATIENTLISTS!
                    else:
                        writer2.writerow(newRow)