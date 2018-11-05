import csv
import random
from copy import deepcopy
from src.HIV.constants import CLADES, REGIONS, AMINOACIDS, POS_ALL
from os.path import join, basename
from os import listdir

from helpers import log_convert, calc_year

AMINO_ACID = 0
PERCENTAGE = 1
DYNAMIC_DATA_FOLDER_NAME = join('data', 'dynamic')
STATIC_DATA_FOLDER_NAME = join('data', 'static')
CONTEMP_PREDICTION_DATA_FOLDER_NAME = join('data', 'contemporary_prediction')
LOG_ZERO_DEFAULT = 0.1


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
            temp = getattr(profs[i], prop_type)
            setattr(profs[i], prop_type, getattr(profs[j], prop_type))
            setattr(profs[j], prop_type, temp)


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
            converted.distr[aa] = log_convert(self.distr[aa])
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
    def __init__(self, amino_acid, clade, region, distr, n_isolates, years, position):
        super().__init__(clade, region, position)
        self.amino_acid = amino_acid  # a single amino acid
        self.years = years
        self.distr = distr  # percentages, in same order as years
        self.numIso = n_isolates  # #isolates, in same order as years
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
            year = calc_year(year)
        for y, distr in zip(self.years, self.distr):
            if y == year:
                return distr
        raise Exception('something is wrong')


# get clade, country, position and return as according enums
def parse_file_name(fn):
    fn = basename(fn)
    # fileName = fileName[fileName.rfind('\\') + 1:]
    [clade, region, position] = fn.split('_')
    position = int(position[:position.rfind('.')])  # remove file extension
    return clade, region, position


# read a file for a clade-region combination into profile instances
def read_dynamic(fn):
    all_profiles = []
    clade, region, position = parse_file_name(fn)

    with open(fn) as file:
        reader = csv.reader(file)

        # read in years
        first_row = next(reader)
        years = []
        for i in range(1, len(first_row)):
            years.append(calc_year(first_row[i]))

        # read in #isolates
        second_row = next(reader)
        n_isolates = []
        for i in range(1, len(second_row)):
            n_isolates.append(int(second_row[i]))

        # read in remaining rows
        for row in reader:
            aa = row[0]
            distr = []
            for i in range(1, len(row)):
                distr.append(float(row[i]))
            profile = DynamicProfile(aa, clade, region, distr, deepcopy(n_isolates), deepcopy(years), position)
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
