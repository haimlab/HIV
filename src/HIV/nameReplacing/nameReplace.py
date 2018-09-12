import csv

def getSpecieKey(s):
    i = s.find('|')
    j = s.find('.', i + 1)
    return s[i + 1:j]

# inputs
matrixFile = 'm.csv'
nameKeyFile = 'n.csv'

# read in the key file
nameKeys = {}
with open(nameKeyFile, 'r') as file:
    for line in file:
        nameKeys[line.strip().split('.')[3]] = line.strip()

# read matrix file and overwrite
with open(matrixFile, 'r') as inFile, open('out.csv', 'w') as outFile:
    reader = csv.reader(inFile)
    writer = csv.writer(outFile, lineterminator = '\n')
    isFirstLine = True
    for line in reader:
        if isFirstLine:
            writer.writerow(line)
            isFirstLine = False
            continue
        try:
            key0 = getSpecieKey(line[0])
            key1 = getSpecieKey(line[1])
            rowToWrite = [nameKeys[key0], nameKeys[key1], line[2]]
        except KeyError:
            rowToWrite = line
        writer.writerow(rowToWrite)