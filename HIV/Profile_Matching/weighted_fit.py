import csv
import numpy as np
import pylab
import sys
import os
from copy import deepcopy
from scipy.optimize import curve_fit
from math import log10
from math import exp
from constants import Clade
from constants import Region
from constants import AminoAcid
from constants import FilterProperties


# object to hold a fit result
class FitResult:
    def __init__(self, slope, y_intercept, r):
        self.y_intercept = y_intercept
        self.slope = slope
        self.r = r

    # calculate year based on percent
    def calcYear(self, percent, default_year):
        try:
            return (percent - self.y_intercept) / self.slope
        except ZeroDivisionError:
            if percent == self.y_intercept:
                return default_year
            else:
                return float('inf')

# object to hold a single distribution profile
class Profile:
    def __init__(self, aminoAcid, clade, region, distr, numIso, years, position):
        self.aminoAcid = aminoAcid # a single amino acid
        self.clade = clade
        self.region = region
        self.years = years
        self.distr = distr # percentages, in same order as years
        self.numIso = numIso # #isolates, in same order as years
        self.position = position
        self.fit = None # a fit object
        self.mostSimilar = None # another profile with minimal euc dist

    # renove data points that have 0 isolates
    def remove_0_isolates(self):
        index = 0
        while index < len(self.numIso):
            if self.numIso[index] == 0:
                del self.numIso[index]
                del self.years[index]
                del self.distr[index]
                continue
            index += 1

    # generate a string to identify this profile
    def tag(self):
        components = [self.clade.value, self.region.value, self.aminoAcid.value]
        return "_".join(components)


# object to hold all profiles and a filtering method
class AllProfiles:
    def __init__(self):
        self.profiles = []

    def add_profile(self, profile):
        self.profiles.append(profile)

    def __filterBy(self, value, property):
        filteredProfiles = AllProfiles()
        for p in self.profiles:
            if property == FilterProperties.AMINOACID:
                curValue = p.aminoAcid
            elif property == FilterProperties.CLADE:
                curValue = p.clade
            elif property == FilterProperties.POSITION:
                curValue = p.position
            elif property == FilterProperties.REGION:
                curValue = p.region
            else:
                raise Exception('Unidentified filter property')
            if curValue == value:
                filteredProfiles.add_profile(p)
        return filteredProfiles


    # filter the profiles according to given criteria
    # returns a result also as AllProfiles instance, so chained filtering can be applied
    def filter(self, clade=None, region=None, aminoAcid=None, position=None):

        # do the filtering
        filtered = self
        if clade is not None:
            filtered = filtered.__filterBy(clade, FilterProperties.CLADE)
        if region is not None:
            filtered = filtered.__filterBy(region, FilterProperties.REGION)
        if aminoAcid is not None:
            filtered = filtered.__filterBy(aminoAcid, FilterProperties.AMINOACID)
        if position is not None:
            filtered = filtered.__filterBy(position, FilterProperties.POSITION)

        return filtered


# calculate r square value of a linear regression, referencing following site
# https://en.wikipedia.org/wiki/Coefficient_of_determination#Definitions
def calc_r_squared(slope, y_intercept, x_data, y_data):

    if len(x_data) != len(y_data):
        raise Exception("input data size mis-match")

    def average(arr):
        total = 0
        for i in arr:
            total += i
        return total / len(arr)

    def sum_squares(arr):
        avg = average(arr)
        var = 0
        for i in arr:
            var += (i - avg) ** 2
        return var

    def sum_residual(slope, y_intercept, x_data, y_data):
        total = 0
        for x, y in zip(x_data, y_data):
            f_i = x * slope + y_intercept
            total += (y - f_i) ** 2
        return total

    sum_squares_tot = sum_squares(y_data)
    sum_squares_residual = sum_residual(slope, y_intercept, x_data, y_data)
    return 1 - sum_squares_residual / sum_squares_tot


# take a profile and compute its weighted linear fit
def calcFit(profile):

    def checkEqual(arr):
        for i in range(0, len(arr) - 1):
            if arr[i] != arr[i + 1]:
                return False
        return True

    profile.remove_0_isolates()

    # to avoid float number round off errors, manually check if all data points are same
    # and assign slope = 0, y_intercept = any data point value, and r square = 1 (perfect fit)
    if checkEqual(profile.distr):
        params = [0, profile.distr[0]]
        r_squared = 1
    else:
        # calculate fit parameters
        sigmas = np.array([1 / n ** .5 for n in profile.numIso]) # weights
        params, cov = curve_fit(lambda x, a, b: a * x + b, profile.years,
                                profile.distr, sigma=sigmas, absolute_sigma=False)
        r_squared = calc_r_squared(params[0], params[1], profile.years, profile.distr)
    profile.fit = FitResult(params[0], params[1], r_squared)


# calculate year as median of the range
# assumes input string to look like "[year1, year2]" with year1 < year2
def calcYear(yearRange):
    commaInd = yearRange.find(',')
    year1 = int(yearRange[1:commaInd])
    year2 = int(yearRange[commaInd + 2:-1])
    return (year1 + year2) / 2


# get clade, country, position and return as according enums
def parseFileName(fileName):
    fileName = fileName[fileName.rfind('\\') + 1:]
    [clade, region, position] = fileName.split('_')
    position = int(position[:position.rfind('.')]) # remove file extension
    return Clade(clade), Region(region), position


# read a file for a clade-region combination into profile instances
def read(fileName):

    allProfiles = []
    clade, region, position = parseFileName(fileName)

    with open(fileName) as file:
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
            profile = Profile(aminoAcid, clade, region, distr, deepcopy(numIso), deepcopy(years), position)
            allProfiles.append(profile)

    return allProfiles


def get_all_profiles(folderName):
    rawFileNames = os.listdir(folderName)
    fileNames = [os.path.join(folderName, fileName) for fileName in rawFileNames]
    profiles = []
    for fileName in fileNames:
        profiles += read(fileName)
    for p in profiles:
        calcFit(p)
    allProfiles = AllProfiles()
    allProfiles.profiles = profiles
    return allProfiles


# calculate euclidean distance, Prof. Haim's approach, p1 and p2 are profiles
def euclideanDist(p1, p2):

    L = 100
    k = -1.2
    x_0 = 5
    logiFunc = lambda x: L / (1 + exp(k * (x - x_0)))  # logistic function
    logConvert = lambda x: 0 if x == 0 else log10(x) + 1  # log convert
    combinedFunc = lambda x: logiFunc(logConvert(x))  # combined logistic and log
    
    # apply transform to data, then square them
    p1 = [combinedFunc(i) for i in p1.distr] 
    p2 = [combinedFunc(i) for i in p2.distr]
    
    # compute euclidean distances
    eucDist = lambda x, y: sum([(a - b) ** 2 for a, b in zip(p1, p2)]) ** .5
    return eucDist(p1, p2)


# # calculate the most similar profiles measured by euclidean distance
# def findMostSimilarProfiles(profiles):
#     for p1 in profiles:
#         if not p1.mostSimilar is None: # skip if similar profiels are already found
#             continue
#         dist = sys.float_info.max
#         for p2 in profiles:
#             if not p1 is p2: # avoid comparing with self
#                 newDist = euclideanDist(p1, p2)
#                 if newDist < dist:
#                     dist = newDist
#                     minDistProfile = p2
#         p1.mostSimilar = minDistProfile
#         minDistProfile.mostSimilar = p1
#
#

def main():
    all_profiles = get_all_profiles(sys.argv[1])
    p = all_profiles.profiles[2]
    x = p.years
    yexact = p.distr
    y = []
    for i in x:
        y.append(p.fit.slope * i + p.fit.y_intercept)
    pylab.plot(x, yexact, 'o', label='Exact')
    pylab.plot(x, y, label='weighted fit')
    pylab.legend(loc='upper right')
    pylab.show()


# actual main
if __name__ == '__main__':
    main()