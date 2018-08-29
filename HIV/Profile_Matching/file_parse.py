import csv
from copy import deepcopy
from constants import FilterProperties, Clade, Region, AminoAcid
from os.path import join
from os import listdir

AMINO_ACID = 0
PERCENTAGE = 1
DATA_FOLDER_NAME = 'data'


class AllStaticProfiles:
    def __init(self):
        self.__profiles = {} # clade -> profile -> region
        # TODO initialize dictionary entries


class StaticProfile:
    def __init__(self, clade, position, region):
        if not isinstance(clade, Clade):
            clade = Clade(clade)
        if not isinstance(region, Region):
            region = Region(region)
        self.__clade = clade
        self.__position = position
        self.__region = region
        self.__distribution = {}

    def clade(self):
        return self.clade

    def position(self):
        return self.position

    def region(self):
        return self.region

    def add_dist(self, amino_acid, percent):
        if not isinstance(amino_acid, AminoAcid):
            amino_acid = AminoAcid(amino_acid)
        self.__distribution[amino_acid] = percent

    def get_distr(self, amino_acid):
        if not isinstance(amino_acid, AminoAcid):
            amino_acid = AminoAcid(amino_acid)
        return self.__distribution[amino_acid]


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


class DynamicProfile:
    def __init__(self, aminoAcid, clade, region, distr, numIso, years, position):
        self.aminoAcid = aminoAcid  # a single amino acid
        self.clade = clade
        self.region = region
        self.years = years
        self.distr = distr  # percentages, in same order as years
        self.numIso = numIso  # #isolates, in same order as years
        self.position = position
        self.fit = None  # a fit object
        self.mostSimilar = None  # another profile with minimal euc dist

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


# get clade, country, position and return as according enums
def parse_file_name(fileName):
    fileName = fileName[fileName.rfind('\\') + 1:]
    [clade, region, position] = fileName.split('_')
    position = int(position[:position.rfind('.')])  # remove file extension
    return Clade(clade), Region(region), position


# read a file for a clade-region combination into profile instances
def read(fileName):

    allProfiles = []
    clade, region, position = parse_file_name(fileName)

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
            profile = DynamicProfile(aminoAcid, clade, region, distr, deepcopy(numIso), deepcopy(years), position)
            allProfiles.append(profile)

    return allProfiles


# read profiles stored in files
def get_all_dynamic_profiles():
    fns = [join(DATA_FOLDER_NAME, fn) for fn in listdir(DATA_FOLDER_NAME)]
    profiles = []
    for fileName in fns:
        profiles += read(fileName)
    allProfiles = AllProfiles()
    allProfiles.profiles = profiles
    return allProfiles


def get_all_static_profiles():
    fns = [join(DATA_FOLDER_NAME, fn) for fn in listdir(DATA_FOLDER_NAME)]
    profiles = []


def get_single_static_profile(file_name):
    clade, region, position = parse_file_name(file_name)
    profile = StaticProfile(clade, position, region)
    with open(file_name) as file:
        reader = csv.reader(file)
        for row in reader:
            profile.add_dist(row[AMINO_ACID], row[PERCENTAGE])
    return profile


# calculate year as median of the range
# assumes input string to look like "[year1, year2]" with year1 < year2
def calcYear(yearRange):
    commaInd = yearRange.find(',')
    year1 = int(yearRange[1:commaInd])
    year2 = int(yearRange[commaInd + 2:-1])
    return (year1 + year2) / 2
