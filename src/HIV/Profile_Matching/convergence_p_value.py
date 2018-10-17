"""
Approach 1 (START by using B_NA)
STEP 1: Calculate the distribution in groups of 100 Envs randomly selected from the years 2007-2015 (repeat 1000 times).
STEP 2: Define centroid of P6 (use all samples from P6).
        STEP 3:  Calculate the distance between each randomly sampled distribution and the centroid of 2007-2015 (altoge
        ther 1,000 distances will be calculated)
STEP 4: Calculate the distance between the historically first 100 Envs (1979-~1992) and P6 centroid.
STEP 5: Calculate ratio:  (number of random P6 distributions for which the distance is greater than P1-P6 distance)/1000
STEP 6: Repeat the above for the second set of 100 Envs (~1993-1998).
"""

from constants import AminoAcid
from csv import reader
from random import sample
from profile_p_value import euc_dist
from argparse import ArgumentParser
from os.path import join
from sys import float_info


FILE_DIR = 'data/convergence'


# compute profile from all given envelopes
# i -> column index of Amino Acids to be counted
def envelopes_to_profile(envs, i):
    count = {aa: 0 for aa in AminoAcid}
    for e in envs:
        count[AminoAcid(e[i])] += 1
    s = sum([count[aa] for aa in count])
    for aa in count:
        count[aa] /= s
    return count


# step 1:
# for n times of shuffle, compute profile of 100 randomly selected samples, collected them
def get_shuffle_profiles(n, fn, pos):
    with open(fn) as f:
        r = reader(f)
        fr = next(r)  # first row
        rows = list(filter(lambda row: 2007 <= int(row[fr.index('Year')]) <= 2015, r))
        return [envelopes_to_profile(sample(rows, 100), fr.index(str(pos))) for _ in range(n)]


# read the nth hundred samples in given position, in historical order, from given envelope table
# and calculate their centroid (which is basically the profile)
def read_nth_group(f_name, pos, n, group_size):
    with open(f_name) as f:
        r = reader(f)
        first_row = next(r)
        for _ in range(0, n):
            for _ in range(0, group_size):
                next(r)
        envs = [next(r) for _ in range(group_size)]
    return envelopes_to_profile(envs, first_row.index(str(pos)))


# get centroid of 2007 to 2015
def get_2007_2015_centroid(fn, pos):
    with open(fn) as f:
        r = reader(f)
        fr = next(r)  # first row
        rows = list(filter(lambda row: 2007 <= int(row[fr.index('Year')]) <= 2015, r))
    return envelopes_to_profile(rows, fr.index(str(pos)))


def main():

    parser = ArgumentParser()
    parser.add_argument('-f', dest='file_name', type=str, required=True)
    parser.add_argument('-g', dest='group_size', type=int, required=True)
    cmd_args = parser.parse_args()
    file_name = join(FILE_DIR, cmd_args.file_name)

    positions = [295, 332, 339, 392, 448]
    for pos in positions:
        for n in range(0, int(float_info.max)):
            rand_prof = get_shuffle_profiles(1000, file_name, pos)  # step 1, get 1000 random shuffled profiles
            cent_07_15 = get_2007_2015_centroid(file_name, pos)  # step 2, calculate centroid of p6
            #  step 3, calculate distance between the 1000 random profiles and 07-15 centroid
            distances = [euc_dist(cent_07_15, i) for i in rand_prof]
            # step 4, distacne between first 100 historical envs with 07-15 cnetroid
            try:
                dist_first_100 = euc_dist(read_nth_group(file_name, pos, n, cmd_args.group_size), cent_07_15)
            except StopIteration:
                break
            # step 5, calculate the ratio as #distance which a random profile from step 1 to 07-15
            # centroid is larger than that of between first_100 envelopes and 07-15 centroid
            ratio = len(list(filter(lambda x: x > dist_first_100, distances))) / 1000
            print(f'position: {pos}, first {n}th p-value: {ratio}')


if __name__ == '__main__':
    main()
