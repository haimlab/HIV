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
        return [envelopes_to_profile(sample(rows, 100), fr.index(str(pos))) for _ in n]


# read the first hundred samples in given position, in historical order, from given envelope table
# and calculate their centroid (which is basically the profile)
def read_first_hundred(f_name, pos):
    with open(f_name) as f:
        r = reader(f)
        first_row = next(r)
        envs = [next(r) for _ in range(100)]
    return envelopes_to_profile(envs, first_row.index(str(pos)))


# get centroid of 2007 to 2015
def get_2007_2015_centroid(fn, pos):
    with open(fn) as f:
        r = reader(f)
        fr = next(r)  # first row
        rows = list(filter(lambda row: 2007 <= int(row[fr.index('Year')]) <= 2015, r))
    return envelopes_to_profile(rows, fr.index(str(pos)))


def main():

    pass

if __name__ == '__main__':
    main()

