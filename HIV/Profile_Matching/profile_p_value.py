from file_parse import get_all_profiles
from math import log10
import constants

LOG_ZERO_DEFAULT = 0.01


def log_transform_all(all_profiles):

    for p in all_profiles.profile:
        p.distr = [logConvert(d) for d in p.distr]


def logConvert(val):
    if val == 0:
        return 0
    else:
        return log10(val) - log10(LOG_ZERO_DEFAULT)


# calculate centroid using mean of each dimension
# output type: {clade: {amino_acid: mean percentage}}
def calc_centroid(all_profiles):
    for clade in constants.Clade:
        for amino_acid in constants.AminoAcid:
            sub_profiles = all_profiles.filter_by(clade=clade, amino_acid=amino_acid)
            for p in sub_profiles.profiles:
                pass


# return {clade, region -> distance}
def quantify_distance(centroids, all_profiles):

    # calculate euclidean distance
    def euclidean_dist(p1, p2):

        # apply transform to data, then square them
        p1 = [logConvert(i) for i in p1.distr]
        p2 = [logConvert(i) for i in p2.distr]

        # compute euclidean distances
        eucDist = lambda x, y: sum([(a - b) ** 2 for a, b in zip(p1, p2)]) ** .5
        return eucDist(p1, p2)


    distances = {}
    for clade in constants.Clade:
        for region in constants.Region:
            sub_profiles = all_profiles.filter_by(clade=clade, region=region)
            if len(sub_profiles.profiles) > 0:  # if the region, clade combo exists
                p = {}


def main():

    # get starting data
    all_profiles = get_all_profiles()
    log_transform_all(all_profiles)

    # grab centroid of each clade

