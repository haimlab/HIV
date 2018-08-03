import csv

# constants
NUM_HEADER_COLS = 4
COUNTRY = 0
YEAR = 1
START = 0
END = 1
POS_WANTED = [295, 332, 339, 392, 448, 662, 663, 664, 665, 667] # epitope positions only
#POS_WANTED = [i for i in range(1, 857)] # all positions
ASIA = ('CN', 'CHINA', 'CN CAT', 'HK', 'IN', 'JP', 'KR', 'MM', 'PH', 'TH', 'TH CAT', 'TW', 'SG')
EU = ('CH', 'CY', 'DE', 'DK', 'ES', 'FR', 'GB', 'IT', 'PL', 'SE', 'NL', 'BE')
NA = ('CA', 'CARR CAT', 'CU', 'DO', 'HT', 'JM', 'TT', 'US', 'US CAT', 'US ICUW')
SA = ('ZA', 'BW') # south africa
ECA = ('DJ', 'BI', 'ET', 'UG', 'KE', 'TZ') # east and central africa, TZ might be removed in future

# inputs
clade = 'C'
# periods = [(1979, 1986), (1987, 1994), (1995, 1999),(2000, 2004), (2005, 2009), (2010, 2015), (2007, 2015)] # both ends inclusive
# countrySet1 = []
# inFileName = '..\\Inputs\\clade_C_all_aligned.csv'
# outDir = '..\\Temp\\'
# end of inputs

# data structure to hold separated files
class Group:
    def __init__(self, periods, countries):
        self.periods = {p: [] for p in periods}
        self.countries = countries
        
    def toStrList(self): # return a list of strings for "pasting" LANL
        strPeriods = []
        for p in self.periods:
            curFileContent = ""
            for row in self.periods[p]:
                curFileContent += row[NUM_HEADER_COLS - 1] + " " # append the identifier
                for i in range(NUM_HEADER_COLS, len(row)):
                    curFileContent += row[i] + " "
                curFileContent += "\n"
            strPeriods.append(curFileContent)
        return strPeriods
                
                
def parseFirstRow(row, pos_wanted):
    pos_inds_wanted = [] # record indices of positions that we want
    for i in range(0, len(row)):
        try:
            pos = int(row[i])
        except ValueError:
            continue
        if pos in pos_wanted:
            pos_inds_wanted.append(i)
    header = []  # create header row
    for i in range(0, len(row)):
        if i < NUM_HEADER_COLS or i in pos_inds_wanted:
            header.append(row[i])
    return pos_inds_wanted, header

def classify(pos_wanted, countryWanted, inFileName, periods):

    # sorted envelopes
    group = Group(periods, countryWanted)

    # classify the file
    with open(inFileName, 'r') as inFile:
        reader = csv.reader(inFile)
        isFirstRow = True
        for row in reader:
            if isFirstRow:
                pos_inds_wanted, header = parseFirstRow(row, pos_wanted)
                isFirstRow = False
                continue
            if row[COUNTRY] in group.countries:  # determine region of envelope
                for period in group.periods:  # determine period of envelope
                    if period[START] <= int(row[YEAR]) <= period[END]:
                        truncatedRow = [row[i] for i in range(0, NUM_HEADER_COLS)]
                        truncatedRow += [row[i] for i in pos_inds_wanted]
                        group.periods[period].append(truncatedRow)
    return group