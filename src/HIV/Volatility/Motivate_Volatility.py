import tkinter
from tkinter import filedialog
import os
from ComputeStandardDevs import getVal
from Bio import SeqIO
import math
from numpy import linalg as LA
import numpy as np

root = tkinter.Tk()
root.withdraw()  # use to hide tkinter window

currdir = os.getcwd()

#names = 'Pat20d1,Pat20d67,Pat1d1,Pat1d110,Pat1d252,Pat2d1,Pat2d156,Pat6d1,Pat6d58,Pat7d1,Pat7d60,Pat9d1,Pat9d160,Pat12d1,Pat12d86,Pat12d114,Pat12d226,Pat13d1,Pat13d113,Pat15d1,Pat15d118,Pat19d1,Pat19d86,Pat3d1,Pat3d172,Pat4d1,Pat4d57,Pat8d1,Pat8d169,Pat10d1,Pat10d57,Pat11d1,Pat11d49,Pat14d1,Pat14d113,Pat16d1,Pat16d87,Pat16d353,Pat17d1,Pat17d41,Pat17d117,Pat18d1,Pat18d107'
#names = 'Pat10d1,Pat11d1,Pat12d1,Pat13d1,Pat14d1,Pat15d1,Pat16d1,Pat17d1,Pat18d1,Pat19d1,Pat1d1,Pat20d1,Pat2d1,Pat3d1,Pat4d1,Pat6d1,Pat7d1,Pat8d1,Pat9d1'
#names = 'Pat20d1,Pat20d67,Pat1d1,Pat1d110,Pat1d252,Pat2d1,Pat2d156,Pat6d1,Pat6d58,Pat7d1,Pat7d60,Pat9d1,Pat9d160,Pat12d1,Pat12d86,Pat12d114,Pat12d226,Pat13d1,Pat13d113,Pat15d1,Pat15d118,Pat19d1,Pat19d86,Pat3d1,Pat3d172,Pat4d1,Pat4d57,Pat8d1,Pat8d169,Pat10d1,Pat10d57,Pat11d1,Pat11d49,Pat14d1,Pat14d113,Pat16d1,Pat16d87,Pat16d353,Pat17d1,Pat17d41,Pat17d117,Pat18d1,Pat18d107'
#names = '701-pre,701-w0,701-w6,701-w7,702-pre,702-w0,702-w5,702-w6,703-pre,703-w5,703-w6,704-pre,704-w0,704-w5,704-w6,707-pre,707-w0,707-w9,707-w10,708-pre,708-w0,708-w9,708-w10,709-pre,709-w5,709-w6,710-pre,710-w19,711-pre,711-w3,711-w4'
names = '701pre,701w0,701w6,701w7,702pre,702w0,702w5,702w6,703pre,703w5,703w6,704pre,704w0,704w5,704w6,707pre,707w0,707w9,707w10,708pre,708w0,708w9,708w10,709pre,709w5,709w6,710pre,710w19,711pre,711w3,711w4'
name = names.split(',')
for nam in name:
    try:
        distances = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Volatility_Data\\V3 Loop Volatility\\For Volatility- Myside\\Distances\\B_700s\\" + nam + '.fas')
        #distances = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Khan_Work\\Maraviroc\\Distances\\" + nam + '.fas')
    except:
        print(nam)
        continue
    # HyPhy Philips Style Matrix has number of Sequences listed in the first position of output#
    # This code brings the matrix into python#
    count = 0
    tot = {''}
    tot.remove('')
    dm = {}

    for row in distances:
        if count == 0:
            count += 1
            continue
        t = row.split()

        fixer = ''
        for i in t[0]:
            if i == '-':
                pass
            else:
                fixer += i
        z0=fixer

        fixer = ''
        for i in t[1]:
            if i == '-':
                pass
            else:
                fixer += i
        z1 = fixer

        if z0 not in tot:
            tot.add(z0)
        if z1 not in tot:
            tot.add(z1)
        dm[(z0, z1)] = t[2]
        dm[(z1, z0)] = t[2]
        count += 1
    numRows = (len(tot))
    if numRows < 2:
        continue
    #################################################################################

    # table = open("F:\\Results\\Vol_B_Eur.xlsx - 2G12.tsv")
    table = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Volatility_Data\\V3 Loop Volatility\\For Volatility- Myside\\Mav_B_700s.txt")
    weights = np.array([1, 1, 1, 1, 1])
    # weights = np.array([1.237176432,0.6845558143,0.2308496402,0.8544987065,0.7873302248]) #2F5
    # weights = np.array([1.981905578,2.017661388,0.712722561,2.232499063,0.5831526609]) #2G12



    # Brings table into Python#
    lines = [[]] * (numRows + 1)
    count = 0

    years = []
    phases = []
    for line in table:
        if len(line) == 0:
            break
        t = line.split()
        if len(t) > 0 and t[0] == 'Patient':
            lines[count] = t
            count += 1
            length = len(t) - 2
            continue
        if len(t) > 1:
            name = 'B' + '.' + t[0] + '.' + t[1]
            if name not in tot:
            #print(name)
                continue
            lines[count] = t
            count += 1
    Positions = [0] * (length)

    c = 0
    for pos in lines[0]:
        if len(pos) != 3:
            continue
        Positions[c] = pos
        c += 1

    count = 0

    num = 0
    vols = [0] * len(Positions)
    cont = False
    posCount = 0
    for pos in Positions:
        EnvIDs = [0] * numRows
        PatNum = [0] * numRows
        positionVals = [0] * (numRows)
        count = 0
        for line in lines:
            if len(line) == 0:
                break
            if count < 1:
                count += 1
                continue
            EnvIDs[count - 1] = 'B' '.' + line[0] + '.' + line[1]
            positionVals[count - 1] = getVal(line[posCount + 2]) / 100 #Depends on the meta data . 2 here for patient and accession
            count += 1
        deltaVals = [[0.0 for x in range(numRows)] for y in range(numRows)]
        volVals = [[0.0 for x in range(numRows)] for y in range(numRows)]
        totalVol = 0
        numCalcs = 1
        for i in range(0, numRows):
            for j in range(0, numRows):
                if i == j:
                    break
                if positionVals[i] == -1000 or positionVals[j] == -1000:
                    continue
                weighted = (positionVals[i] - positionVals[j])
                deltaVals[i][j] = math.pow(weighted, 2)
                try:
                    t = float(dm[(EnvIDs[i], EnvIDs[j])])
                    r = round(float(dm[(EnvIDs[i], EnvIDs[j])]), 8)
                    if r > 0:
                        volVals[i][j] = round(deltaVals[i][j] / float(dm[(EnvIDs[i], EnvIDs[j])]), 8)
                        totalVol += volVals[i][j]
                        numCalcs += 1
                except:
                    pass
        vols[posCount] = str(totalVol / numCalcs)
        posCount += 1
        num = count - 1
    # if cont:
    #     continue
    vol_string = " "
    for vol in vols:
        vol_string += str(vol) + ' '
    try:
        print(str(num) + ' ' + line[0] + ' ' + line[1] + vol_string)
    except:
        print('FAILURE')

