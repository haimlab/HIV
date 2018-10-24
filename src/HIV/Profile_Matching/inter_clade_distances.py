from argparse import ArgumentParser
from csv import reader
from convergence_p_value import envelopes_to_profile
from profile_p_value import euc_dist
from constants import POS_2G12


# take a 1D list of years and put them into groups of two
def format_period_ranges(years):
    if len(years) % 2 != 0:
        raise Exception('cannot group odd length lists')
    return [(years[i], years[i + 1]) for i in range(len(years) - 1)]


# filter an input file based on year
def filter_by_year(low, high, file_name):
    with open(file_name) as f:
        r = reader(f)
        year_ind = next(r).index('Year')
        return list(filter(lambda r: low <= r[year_ind] <= high, [row for row in r]))


# assume
# - each file has its own list of year_ranges
# - all such lists have the same length
def main():
    file_names = []
    year_ranges = format_period_ranges([
        1987, 1999,
        2000, 2002,
        2003, 2004,
        2005, 2005,
        2006, 2007,
        2008, 2009,
        2010, 2015
    ])
    centroids = {}

    # compute all centroids in each period
    for p in POS_2G12:
        for fn, (low, high) in zip(file_names, year_ranges):
            with open(fn) as f:
                aa_ind = next(reader(f)).index(p)
            env = envelopes_to_profile(filter_by_year(low, high, fn), aa_ind)
            centroids[(low, high)][fn] = env

    # compute distances between all centroids of the same period
    distances = {y_range: {} for y_range in centroids}
    for y_range in centroids:
        for f_1 in centroids[y_range]:
            for f_2 in centroids[y_range]:
                key = min(f_1, f_2), min(f_1, f_2)
                distances[y_range][key] = euc_dist(centroids[y_range][f_1], centroids[y_range][f_2])


    # print all the results, for now
    for y_range in distances:
        print(f'year_range: {y_range}')
        for key_1, key_2 in y_range:
            print(f'distance between {key_1} and {key_2} is')
