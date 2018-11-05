"""
This file has code for calculating clade and position specificity.

Clade Specificity Steps:
1. Cent = Location of position centroid (calculated using all regions in all clades for that position).
2. D(Cent , Reg&Clade) = distance between profile of the centroid of each position and the profile of each region&clade.
3. D(Cent , Cent): Distance between the centroids of all positions.
4. R = [Average D(Cent , Reg&Clade)] DIVIDED BY Average_D(cent, cent)
5. Shuffle the position numbers for all region/clade profiles and calculate the above again.
6. The p value is the number of times (out of 10,000 repeats) that the ratio R of a given position is smaller for the
   shuffled than for the results of the test using the correctly-labeled data.

Posiitonal Specificity Steps:
Similar to those of clade specificity, but position labels are shuffled, and subgroups are grouped by positions.
"""

from src.HIV.Profile_Matching.file_parse import get_all_static_profiles, AllStaticProfiles
from numpy import asarray
from scipy.cluster.vq import kmeans
from src.HIV.constants import AMINOACIDS, POS_PNGS, POS_2F5, POS_2G12
from argparse import ArgumentParser
from copy import deepcopy
from src.HIV.Profile_Matching.helpers import euc_dist


def calc_all_centroids(all_profiles, prop_type):
    centroids = {}
    for p in all_profiles.attr_list(prop_type):
        sub = all_profiles.filter(p)
        centroid = calc_centroid(sub)
        centroids[p] = centroid
    return centroids  # {prop_type - {amino acid -> percentage}


# calculate centroid by using only 1 centroid for k-means clustering
# output type: {clade: {amino_acid: mean percentage}}
def calc_centroid(all_profiles):
    centroid, _ = kmeans(asarray([[prof.distr[aa] for aa in AMINOACIDS] for prof in all_profiles.profiles]), 1)

    # format output
    centroid = centroid[0]
    centroid_dict = {}
    index = 0
    for aa in AMINOACIDS:
        centroid_dict[aa] = centroid[index]
        index += 1
    return centroid_dict


def clade_specificity_one_round(profs):

    # compute distance within
    centroids = calc_all_centroids(profs, 'clade')
    distances = []
    for prop in centroids:
        sub_profiles = profs.filter(prop)
        c = centroids[prop]
        for p in sub_profiles.profiles:
            distances.append(euc_dist(c, p.distr))
    dist_within = sum(distances) / len(distances)

    # compute distance without
    c_list = [centroids[prop] for prop in centroids]
    distances = []
    for i in range(0, len(c_list)):
        for j in range(i + 1, len(c_list)):
            distances.append(euc_dist(c_list[i], c_list[j]))
    dist_without = sum(distances) / len(distances)

    return dist_within / dist_without


def clade_specificity(num_shuffle, all_profiles, positions):
    for pos in positions:
        sub = deepcopy(all_profiles.filter(pos))
        std_rat = clade_specificity_one_round(sub)
        num_above = 0
        num_below = 0
        num_equal = 0
        for i in range(0, num_shuffle):
            sub.shuffle('clade')
            r = clade_specificity_one_round(sub)
            if r > std_rat:
                num_above += 1
            elif r < std_rat:
                num_below += 1
            else:
                num_equal += 1
        print(f'position: {pos}, p-value: {num_below / num_above}')
        print(f'below: {num_below}, above: {num_above}, equal: {num_equal}')


def select_sub_group(raw, clade_region_pairs, positions):
    pairs = [i for i in (p.split(',') for p in clade_region_pairs)]
    filters = []
    for p in positions:
        filters += [(p, c, r) for [c, r] in pairs]
    all_prof = AllStaticProfiles()
    for f in filters:
        sub = raw.filter(*f).profiles
        for s in sub:
            if s.clade == 'C' and s.position == 295:
                continue
            if s.clade == 'AE' and (s.position == 332 or s.position == 339):
                continue
            all_prof.profiles.append(s)
    return all_prof


# input profiles should already be shuffled
def position_specificity_one_round(all_profiles, cur_prop_val):

    # calculate distance without (distance relative to other sub groups)
    all_centroids = calc_all_centroids(all_profiles, 'position')
    cur_centroid = all_centroids[cur_prop_val]
    del all_centroids[cur_prop_val]
    dist_without = sum([euc_dist(cur_centroid, all_centroids[c]) for c in all_centroids]) / len(all_centroids)

    # calculate distance within
    sub_profiles = all_profiles.filter(cur_prop_val)
    total = sum([euc_dist(cur_centroid, p.distr) for p in sub_profiles.profiles])
    dist_within = total / len(sub_profiles.profiles)

    return dist_within / dist_without


def position_specificity(num_shuffle, all_profiles):

    positions = all_profiles.attr_list('position')
    shuffle_copy = deepcopy(all_profiles)

    for p in positions:
        std_rat = position_specificity_one_round(all_profiles, p)
        num_above = 0
        num_below = 0
        for i in range(0, num_shuffle):
            shuffle_copy.shuffle('position')
            shuffled_rat = position_specificity_one_round(shuffle_copy, p)
            if shuffled_rat > std_rat:
                num_above += 1
            elif shuffled_rat < std_rat:
                num_below += 1
        try:
            print(f'{p}: p-value: {num_below / num_above}')
            print(f'times shuffled: {num_shuffle}, un-shuffled ratio: {std_rat}')
            print(f'above: {num_above}, below: {num_below}, dropped: {num_shuffle - num_above - num_below}')
        except ZeroDivisionError:
            print(f'{p}: p-value: zero division')


def main():
    parser = ArgumentParser()
    parser.add_argument('-t', dest='type', type=str, required=True)
    parser.add_argument('-n', dest='num_shuffle', type=int, required=True)
    parser.add_argument('-e', dest='epitope', type=str, required=True)
    parser.add_argument('-p', nargs='+', dest='clade_region_pairs', required=True)
    cmd_args = parser.parse_args()

    if cmd_args.epitope == '2g12':
        positions = POS_2G12
    elif cmd_args.epitope == '2f5':
        positions = POS_2F5
    elif cmd_args.epitope == 'pngs':
        positions = POS_PNGS
    else:
        raise Exception('invalid epitope type')

    all_prof = select_sub_group(get_all_static_profiles(), cmd_args.clade_region_pairs, positions)
    all_prof = all_prof.log_convert()

    if cmd_args.type == 'clade':
        clade_specificity(cmd_args.num_shuffle, all_prof, positions)
    elif cmd_args.type == 'position':
        position_specificity(cmd_args.num_shuffle, all_prof)
        print(positions)
    else:
        raise Exception('invalid')


if __name__ == '__main__':
    main()
