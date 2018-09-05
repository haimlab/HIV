# Split a large clade fasta file into individual smaller patient-specific ones

import tkinter
from tkinter import filedialog
import os
from Bio import SeqIO

# obtain input/output directory and files
root = tkinter.Tk()
root.withdraw()
curDir = os.getcwd()
tempDir = filedialog.askopenfilename(parent = root, initialdir = curDir, \
                                     title='Feature Table')
fasta_sequences = SeqIO.parse(open(tempDir),'fasta')
outDir = filedialog.askdirectory()

# read all data into memory
allEntries = {}
for fasta in fasta_sequences:
    try:
        curID = fasta.id.replace('/', '_').split('.')[2]
    except IndexError:
        continue
    if curID not in allEntries:
        allEntries[curID] = ''
    allEntries[curID] += ('>' + str(fasta.id) + '\n' + str(fasta.seq) +'\n')

# write all data after finished, and a file with just patient names
allNames = ''
for curID in allEntries:
    allNames += curID + ","
    path = os.path.join(outDir, curID + '.fas')
    file = open(path, 'w')
    file.write(allEntries[curID])
    file.close()
allNames = allNames[:len(allNames) - 1] # remove last comma
path = os.path.join(outDir, 'allPatientNames.txt')
with open(path, 'w') as file:
    file.write(allNames)
