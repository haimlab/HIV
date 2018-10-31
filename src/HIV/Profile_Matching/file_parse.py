import csv
import random
from copy import deepcopy
from constants import CLADES, REGIONS, AMINOACIDS, POS_ALL
from os.path import join, basename
from os import listdir
from math import log10
from functools import reduce

AMINO_ACID = 0
PERCENTAGE = 1
DYNAMIC_DATA_FOLDER_NAME = join('data', 'dynamic')
STATIC_DATA_FOLDER_NAME = join('data', 'static')
CONTEMP_PREDICTION_DATA_FOLDER_NAME = join('data', 'contemporary_prediction')
LOG_ZERO_DEFAULT = 0.1


# TODO add option of parsing selected subset of input files only, given clade, region, position
# TODO add distance member method to profile class
def logConvert(val):
    val = log10(LOG_ZERO_DEFAULT) if val < LOG_ZERO_DEFAULT else log10(val)
    return val - log10(LOG_ZERO_DEFAULT)


class Profile:
    """ abstract generic class for static and dynamic profiles """

    def __init__(self, clade, region, position):
        self.clade = clade
        self.position = position
        self.region = region


class AllProfiles:
    """ generic abstract class for static and dynamic profile managers """

    def __init__(self):
        self.profiles = []

    # filter the profiles according to given criteria
    # returns a result also as AllProfiles instance, so chained filtering can be applied
    def filter(self, *args):
        # verify if filter arguments are valid
        invalid_args = list(filter(lambda x: x not in CLADES + REGIONS + POS_ALL + AMINOACIDS, args))
        if len(invalid_args) > 0:
            invalid_args = list(map(str, invalid_args))
            raise Exception('invalid filters:\n' + '\n'.join(invalid_args))

        filtered = AllProfiles()
        for p in self.profiles:
            val_set = {p.clade, p.region, p.position}
            if hasattr(p, 'amino_acid'):
                val_set.add(p.amino_acid)
            if all([arg in val_set for arg in args]):
                filtered.profiles.append(p)
        return filtered

    def get_only_profile(self):
        if len(self.profiles) != 1:
            raise Exception('cannot find single profile')
        else:
            return self.profiles[0]

    def attr_list(self, prop_type):
        return list({getattr(p, prop_type) for p in self.profiles})

    def shuffle(self, prop_type):
        # so that we don't shuffle things back to previous orders when
        # we do two multiple shuffles in a row
        random.seed()
        profs = self.profiles
        for i in range(len(profs) - 1):
            j = random.randint(i, len(profs) - 1)
            if prop_type == 'clade':
                temp = profs[i].clade
                profs[i].clade = profs[j].clade
                profs[j].clade = temp
            elif prop_type == 'position':
                temp = profs[i].position
                profs[i].position = profs[j].position
                profs[j].position = temp
            elif prop_type == 'region':
                temp = profs[i].region
                profs[i].region = profs[j].region
                profs[j].region = temp
            else:
                raise Exception('Unimplemented shuffle property')


class AllStaticProfiles(AllProfiles):
    def __init__(self):
        super().__init__()

    def log_convert(self):
        converted = AllStaticProfiles()
        for p in self.profiles:
            converted.profiles.append(p.log_convert())
        return converted


class StaticProfile(Profile):
    def __init__(self, clade, region, position):
        super().__init__(clade, region, position)
        self.distr = {}  # amino acid -> percent

    def log_convert(self):
        converted = StaticProfile(self.clade, self.region, self.position)
        for aa in self.distr:
            converted.distr[aa] = logConvert(self.distr[aa])
        return converted


class AllDynamicProfiles(AllProfiles):
    def __init__(self):
        super().__init__()

    def get_profile(self, clade, region, position, year):
        prof = []
        p = self.filter(clade, region, position)
        for aa in AMINOACIDS:
            i = p.filter(aa).get_only_profile()
            prof.append(i.get_distr(year))
        return prof


class DynamicProfile(Profile):
    def __init__(self, amino_acid, clade, region, distr, numIso, years, position):
        super().__init__(clade, region, position)
        self.amino_acid = amino_acid  # a single amino acid
        self.years = years
        self.distr = distr  # percentages, in same order as years
        self.numIso = numIso  # #isolates, in same order as years
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

    def get_distr(self, year):
        if type(year) is str:
            year = calcYear(year)
        for y, distr in zip(self.years, self.distr):
            if y == year:
                return distr
        raise Exception('something is wrong')


# get clade, country, position and return as according enums
def parse_file_name(fileName):
    fileName = basename(fileName)
    # fileName = fileName[fileName.rfind('\\') + 1:]
    [clade, region, position] = fileName.split('_')
    position = int(position[:position.rfind('.')])  # remove file extension
    return clade, region, position


# read a file for a clade-region combination into profile instances
def read_dynamic(fileName):
    all_profiles = []
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
            aa = row[0]
            distr = []
            for i in range(1, len(row)):
                distr.append(float(row[i]))
            profile = DynamicProfile(aa, clade, region, distr, deepcopy(numIso), deepcopy(years), position)
            all_profiles.append(profile)

    return all_profiles


# read profiles stored in files
def get_all_dynamic_profiles():
    fns = [join(DYNAMIC_DATA_FOLDER_NAME, fn) for fn in listdir(DYNAMIC_DATA_FOLDER_NAME)]
    all_profs = AllDynamicProfiles()
    for fileName in fns:
        profs = read_dynamic(fileName)
        for p in profs:
            all_profs.profiles.append(p)
    return all_profs


def get_all_static_profiles():
    fns = [join(STATIC_DATA_FOLDER_NAME, fn) for fn in listdir(STATIC_DATA_FOLDER_NAME)]
    all_profiles = AllStaticProfiles()
    for fn in fns:
        all_profiles.profiles.append(read_static(fn))
    return all_profiles


def get_all_contemporary_prediction_profiles():
    fns = [join(CONTEMP_PREDICTION_DATA_FOLDER_NAME, fn) for fn in listdir(CONTEMP_PREDICTION_DATA_FOLDER_NAME)]
    all_profs = AllDynamicProfiles()
    for fileName in fns:
        profs = read_dynamic(fileName)
        for p in profs:
            all_profs.profiles.append(p)
    return all_profs


def read_static(file_name):
    clade, region, position = parse_file_name(file_name)
    profile = StaticProfile(clade, region, position)
    with open(file_name) as file:
        reader = csv.reader(file)
        for row in reader:
            profile.distr[row[AMINO_ACID]] = float(row[PERCENTAGE])
    return profile


# calculate year as median of the range
# assumes input string to look like "[year1, year2]" with year1 < year2
def calcYear(yearRange):
    commaInd = yearRange.find(',')
    year1 = int(yearRange[1:commaInd])
    year2 = int(yearRange[commaInd + 2:-1])
    return (year1 + year2) / 2


if __name__ == '__main__':
    a = get_all_static_profiles()
