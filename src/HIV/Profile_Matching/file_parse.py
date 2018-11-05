import csv
import random
from src.HIV.constants import CLADES, REGIONS, AMINOACIDS, POS_ALL
from os.path import join, basename
from os import listdir
from src.HIV.Profile_Matching.helpers import log_convert

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


# get clade, country, position and return as according enums
def parse_file_name(fn):
    fn = basename(fn)
    [clade, region, position] = fn.split('_')
    position = int(position[:position.rfind('.')])  # remove file extension
    return clade, region, position


def get_all_static_profiles():
    fns = [join(STATIC_DATA_FOLDER_NAME, fn) for fn in listdir(STATIC_DATA_FOLDER_NAME)]
    all_profiles = AllStaticProfiles()
    for fn in fns:
        all_profiles.profiles.append(read_static(fn))
    return all_profiles


def read_static(file_name):
    clade, region, position = parse_file_name(file_name)
    profile = StaticProfile(clade, region, position)
    with open(file_name) as file:
        reader = csv.reader(file)
        for row in reader:
            profile.distr[row[AMINO_ACID]] = float(row[PERCENTAGE])
    return profile
