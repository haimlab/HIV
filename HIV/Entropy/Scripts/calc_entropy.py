import requests
import region_separation
from html.parser import HTMLParser
import csv

# script level constants
ENTROPY1_CALC_URL = 'https://www.hiv.lanl.gov/cgi-bin/ENTROPY/entropy_one.cgi'
LANL_HOME_URL = 'https://www.hiv.lanl.gov'
NUM_HEADER_COLS = 4
COUNTRY = 0
YEAR = 1
START = 0
END = 1
POS_WANTED = [295, 332, 339, 392, 448, 662, 663, 664, 665, 667]  # epitope positions only
periods = [
    (1979, 1986),  # both ends inclusive
    (1987, 1994),
    (1995, 1999),
    (2000, 2004),
    (2005, 2009),
    (2010, 2015),
    (2007, 2015)
]
ASIA = ('CN', 'CHINA', 'CN CAT', 'HK', 'IN', 'JP', 'KR', 'MM', 'PH', 'TH', 'TH CAT', 'TW', 'SG')
EU = ('CH', 'CY', 'DE', 'DK', 'ES', 'FR', 'GB', 'IT', 'PL', 'SE', 'NL', 'BE')  # europe
NA = ('CA', 'CARR CAT', 'CU', 'DO', 'HT', 'JM', 'TT', 'US', 'US CAT', 'US ICUW')  # north americ
SA = ('ZA', 'BW')  # south africa
CLADE_AE_FILE_NAME = '..\\Inputs\\clade_AE_all_aligned.csv'
CLADE_B_FILE_NAME = '..\\Inputs\\clade_B_all_aligned.csv'
CLADE_C_FILE_NAME = '..\\Inputs\\clade_C_all_aligned.csv'

# inputs
clade = 'C'
outDir = '..\\Temp\\'


# determine which pre-stored input file to use
def choose_input(clade):
    if clade == 'AE':
        return CLADE_AE_FILE_NAME
    if clade == 'B':
        return CLADE_B_FILE_NAME
    if clade == 'C':
        return CLADE_C_FILE_NAME
    else:
        raise Exception('Unknown Clade')


# find the index of the nth occurence of target in string
def findNth(string, target, n):
    if n < 1:
        raise Exception('Error')

    start = -1
    while n > 0:
        start = string.find(target, start + 1)
        n -= 1
    return start


inFileName = choose_input(clade)
group = region_separation.classify(POS_WANTED, SA, inFileName, periods)
strList = group.toStrList()

form_data = {'seq_input': strList[2]}

r = requests.post(ENTROPY1_CALC_URL, data=form_data)  # calcualte enetropy
src = r.text[findNth(r.text, '<', 3) + 12: r.text.find('>', findNth(r.text, '<', 3)) - 1]

r = requests.get(LANL_HOME_URL + src)  # go to entropy page
i = findNth(r.text, '<', 22)
j = r.text.find('>', i)
textFormUrl = r.text[i + 9: j - 20]
r = requests.get(LANL_HOME_URL + textFormUrl)  # go to entropy text form page
htmlText = r.text
with open('sample.html', 'w') as sampleFile:
    sampleFile.write(r.text)


# parse html
class MyHTMLParser(HTMLParser):

    def error(self, message):
        raise Exception("HTML parsing error")

    def __init__(self):
        super().__init__()
        self.colNum = 3
        self.curColNum = 1
        self.table = []
        self.row = []
        self.insideTable = False

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.insideTable = True

    def handle_endtag(self, tag):
        if tag == 'table':
            self.insideTable = False

    def handle_data(self, data):
        if self.insideTable:
            if self.insideTable:
                self.row.append(str(data))
                self.curColNum += 1
                if self.curColNum > self.colNum:
                    self.table.append(self.row)
                    self.row = []
                    self.curColNum = 1
        print("Data     :", data)


parser = MyHTMLParser()
parser.feed(htmlText)

# write parsed table to output file
with open('..\\Temp\\out.csv', 'w') as outFile:
    writer = csv.writer(outFile, lineterminator='\n')
    writer.writerows(parser.table)
