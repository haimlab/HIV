import csv
import numpy as np
from scipy.optimize import curve_fit
import pylab
from copy import deepcopy
from math import log10
from math import exp
from constants import Clade
from constants import Region
from constants import AminoAcid


# object to hold a fit result
class FitResult:
    def __init__(self, slope, y_intercept, r):
        self.y_intercept = y_intercept
        self.slope = slope
        self.r = r


# object to hold a single distribution profile
class Profile:
    def __init__(self, aminoAcid, clade, region, distr, numIso, years):
        self.aminoAcid = aminoAcid
        self.clade = clade
        self.region = region
        self.years = years
        self.distr = distr # percentages, in same order as years
        self.numIso = numIso # #isolates, in same order as years
        self.fit = None # a fit object
        
    # generate a string to identify this profile
    def tag(self):
        components = [self.clade.value, self.region.value, self.aminoAcid.value]
        return "_".join(components)

# take a profile and compute its weighted linear fit
def calcFit(profile): 
    sigmas = np.array([1 / n ** .5 for n in profile.numIso]) # weights
    params, cov = curve_fit(lambda x, a, b: a * x + b, profile.years,
                            profile.distr, sigma=sigmas, absolute_sigma=False)
    profile.fit = FitResult(params[0], params[1], cov)
    
    
# calculate year as median of the range
# assumes input string to look like "[year1, year2]" with year1 < year2
def calcYear(yearRange):
    commaInd = yearRange.find(',')
    year1 = int(yearRange[1:commaInd])
    year2 = int(yearRange[commaInd + 2:-1])
    return (year1 + year2) / 2


# get clade, country, position and return as according enums
def parseFileName(fileName):
    [clade, region, position] = fileName.split('_')
    position = int(position[:position.rfind('.')]) # remove file extension
    return Clade(clade), Region(region), position


# read a file for a clade-region combination into profile instances
def read(fileName):
    
    allProfiles = []
    clade, region, position = parseFileName(fileName)
    
    with open(fileName, 'r') as file:
        reader = csv.reader(file)
        
        # read in years
        firstRow = next(reader)
        years = []
        for i in range(1, len(firstRow)):
            years.append(calcYear(firstRow[i]))
        
        # read in #isolates
        secondRow = next(reader)
        numIso = []
        for i in range(1, len(secondRow)):
            numIso.append(int(secondRow[i]))
        
        # read in remaining rows
        for row in reader:
            aminoAcid = AminoAcid(row[0])
            distr = []
            for i in range(1, len(row)):
                distr.append(float(row[i]))
            profile = Profile(aminoAcid, clade, region, distr, numIso, years)
            allProfiles.append(profile)
        
    return allProfiles


# calculate euclidean distance, Prof. Haim's approach, p1 and p2 are profiles
def euclideanDist(p1, p2):
    
    L = 4
    k = -1.5
    x_0 = 3.5
    logiFunc = lambda x: L / 1 + exp(k * (x - x_0)) # logistic function
    logConvert = lambda x: 0 if x == 0 else log10(x) + 1 # log convert
    combinedFunc = lambda x: logiFunc(logConvert(x)) # combined logistic and log
    
    # apply transform to data, then square them
    p1 = [combinedFunc(i) for i in p1.distr] 
    p2 = [combinedFunc(i) for i in p2.distr]
    
    # compute euclidean distances
    eucDist = lambda x, y: sum([(a - b) ** 2 for a, b in zip(p1, p2)]) ** .5
    return eucDist(p1, p2)
    
    
if __name__ == '__main__':
    
    profiles = read('B_EU_295.csv')
    for p in profiles:
        calcFit(p)
    
    # plot second profile for verification
    p = profiles[0]
    x = p.years
    yexact = p.distr
    y = []
    for i in x:
        y.append(p.fit.slope * i + p.fit.y_intercept)
    pylab.plot(x, yexact, 'o', label='Exact')
    pylab.plot(x, y, label='weighted fit')
    pylab.legend(loc='upper right')
    pylab.show()
    
    # verifiy the transform functions
    L = 4
    k = -1.5
    x_0 = 3.5
    logiFunc = lambda x: L / (1 + exp(k * (x - x_0))) # logistic function
    logConvert = lambda x: 0 if x == 0 else log10(x) + 1 # log convert
    combinedFunc = lambda x: logiFunc(logConvert(x)) # combined logistic and log
    # use B NA / all B
    for p in profiles:
        pylab.figure(p.tag())
        y = [combinedFunc(i) for i in p.distr]
        print(len(y))
        pylab.plot(p.distr, y)
        pylab.plot(p.distr, y, 'o')
        pylab.xlabel("actual")
        pylab.ylabel("converted")
        pylab.title(p.tag())
        