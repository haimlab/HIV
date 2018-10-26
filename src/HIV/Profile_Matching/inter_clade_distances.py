from csv import reader
from helpers import envelopes_to_profile
from profile_p_value import euc_dist
from constants import POS_2G12, AminoAcid
from os.path import join, basename
from random import sample
from scipy.cluster.vq import kmeans
from numpy import asarray


# take a 1D list of years and put them into groups of two
def format_period_ranges(years):
    if len(years) % 2 != 0:
        raise Exception('cannot group odd length lists')
    return [(years[i], years[i + 1]) for i in range(0, len(years) - 1, 2)]


# filter an input file based on year
def filter_by_year(low, high, file_name):
    with open(file_name) as f:
        r = reader(f)
        year_ind = next(r).index('Year')
        return list(filter(lambda r: low <= int(r[year_ind]) <= high, [row for row in r]))


# assume
# - each file has its own list of year_ranges
# - all such lists have the same length
def main():
    file_dir = join('data', 'inter_clade')
    file_names = ['B_NA&EU.csv', 'C_SA&ECA&EU.csv']
    num_samples = [
        (260, 133),
        (114, 144),
        (136, 266),
        (107, 125),
        (210, 128),
        (185, 177),
        (170, 165)
    ]
    year_ranges = format_period_ranges([
        1987, 1999,
        2000, 2002,
        2003, 2004,
        2005, 2005,
        2006, 2007,
        2008, 2009,
        2010, 2015
    ])
    print(year_ranges)
    centroids = {p: {y: {} for y in year_ranges} for p in POS_2G12}

    # compute all centroids in each period
    for p in POS_2G12:
        for (low, high), (n_b, n_c) in zip(year_ranges, num_samples):
            for fn in file_names:
                n_sample = min(n_b, n_c)
                fn = join(file_dir, fn)
                with open(fn) as f:
                    aa_ind = next(reader(f)).index(str(p))
                centroids_rand = []
                candidate_profs = filter_by_year(low, high, fn)
                print(f'position {p}, range: {low} - {high}, file {fn}')
                for i in range(10000):
                    centroids_rand.append(envelopes_to_profile(sample(candidate_profs, n_sample), aa_ind))
                centroids_list_rand = [[c[aa] for aa in AminoAcid] for c in centroids_rand]
                cent_as_list = kmeans(asarray(centroids_list_rand), 1)[0][0]
                cent = {aa: i for aa, i in zip(AminoAcid, cent_as_list)}
                centroids[p][(low, high)][basename(fn)] = cent

    # compute distances between all centroids of the same period
    distances = {p: {y_range: {} for y_range in centroids[p]} for p in POS_2G12}
    for p in POS_2G12:
        for y_range in centroids[p]:
            for f_1 in centroids[p][y_range]:
                for f_2 in centroids[p][y_range]:
                    if f_1 == f_2:
                        continue
                    dist = euc_dist(centroids[p][y_range][f_1], centroids[p][y_range][f_2])
                    f_1 = basename(f_1)
                    f_2 = basename(f_2)
                    key = min(f_1, f_2), max(f_1, f_2)
                    distances[p][y_range][key] = dist

    # print all the results, for now
    for p in POS_2G12:
        print(f'\nposition: {p}')
        for key_1, key_2 in distances[p]:
            print(f'{key_1} - {key_2}')
            for f_1, f_2 in distances[p][(key_1, key_2)]:
                print(f'{f_1} - {f_2}: {distances[p][(key_1, key_2)][(f_1, f_2)]}')


if __name__ == '__main__':
    main()
