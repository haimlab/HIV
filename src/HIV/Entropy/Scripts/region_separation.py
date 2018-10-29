import csv

# constants
NUM_HEADER_COLS = 4
COUNTRY = 0
YEAR = 1
START = 0
END = 1


# data structure to hold separated files
class Group:
    def __init__(self, periods, countries):
        self.periods = {p: [] for p in periods}
        self.countries = countries

    def toStrList(self):  # return a list of strings for "pasting" LANL
        strPeriods = []
        for p in self.periods:
            curFileContent = ""
            for row in self.periods[p]:
                curFileContent += row[NUM_HEADER_COLS - 1] + " "  # append the identifier
                for i in range(NUM_HEADER_COLS, len(row)):
                    curFileContent += row[i] + " "
                curFileContent += "\n"
            strPeriods.append(curFileContent)
        return strPeriods


# record indices of positions wanted with respect to first column of file
def parseFirstRow(row, pos_wanted):
    pos_inds_wanted = []
    for i in range(0, len(row)):
        try:
            pos = int(row[i])
        except ValueError:
            continue
        if pos in pos_wanted:
            pos_inds_wanted.append(i)
    return pos_inds_wanted


def classify(pos_wanted, countryWanted, inFileName, periods):
    group = Group(periods, countryWanted)  # object holding classification results
    with open(inFileName, 'r') as inFile:  # classify the file
        reader = csv.reader(inFile)
        pos_inds_wanted = parseFirstRow(next(reader), pos_wanted)
        for row in reader:
            if row[COUNTRY] in group.countries:  # determine region of envelope
                for period in group.periods:  # determine period of envelope
                    if period[START] <= int(row[YEAR]) <= period[END]:
                        truncatedRow = [row[i] for i in range(0, NUM_HEADER_COLS)]
                        truncatedRow += [row[i] for i in pos_inds_wanted]
                        group.periods[period].append(truncatedRow)
    return group
