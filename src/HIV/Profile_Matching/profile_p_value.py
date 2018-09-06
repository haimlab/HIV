from file_parse import get_all_static_profiles
from numpy import asarray
from scipy.cluster.vq import kmeans
from math import log10
import constants
from constants import FilterProperties


# output type {clade - {amino acid -> percentage}
def calc_all_centroids(all_profiles):
    centroids = {}
    clades = all_profiles.attr_list(FilterProperties.CLADE)
    for clade in clades:
        sub = all_profiles.filter(clade=clade)
        centroid = calc_centroid(sub)
        centroids[clade] = centroid
    return centroids


# calculate centroid using mean of each dimension
# delegates to kmeans to scipy
# output type: {clade: {amino_acid: mean percentage}}
def calc_centroid(all_profiles):

    # format input for scipy
    vectors = []
    dims = None
    for prof in all_profiles.get_all_profiles():
        cur_dim = prof.dim()
        if dims is None:
            dims = cur_dim
        if cur_dim != dims:  # different amino acid in profiles
            raise Exception('profiles contain different amino acids')
        vectors.append([prof.get_distr(aa) for aa in dims])

    # compute centroid using scipy
    centroid, _ = kmeans(asarray(vectors), 1)

    # format output
    centroid = centroid[0]
    centroid_dict = {}
    index = 0
    for aa in dims:
        centroid_dict[aa] = centroid[index]
        index += 1
    return centroid_dict


# calculate euclidean distance
def euc_dist(p1, p2):
    if p1.keys() != p2.keys():
        raise Exception('dimension mismatch')
    p1 = [p1[key] for key in p1.keys()]
    p2 = [p2[key] for key in p2.keys()]
    return (sum([(a - b) ** 2 for a, b in zip(p1, p2)])) ** .5


# return {clade, region -> distance}
# centroids: {clade -> {aa -> centroid_val}}
def quantify_distances(centroids, all_profiles):
    distances = {}
    for clade in centroids:
        sub_profiles = all_profiles.filter(clade=clade)
        centroid = centroids[clade]
        for prof in sub_profiles.get_all_profiles():
            distr = prof.get_entire_distr()
            distances[(clade, prof.region())] = euc_dist(centroid, distr)
    return distances


# average distance between all centroids
# input {clade -> {aa -> centroid_val}}
def avg_centroid(centroids):
    c_list = [centroids[clade] for clade in centroids]
    dist_list = []
    for i in range(0, len(c_list)):
        c = c_list[i]
        for j in range(0, len(c_list)):
            if j != i:
                dist_list.append(euc_dist(c, c_list[j]))
    return sum(dist_list) / len(dist_list)


# calculates the ratio of average distance between all centroids and average distance between all centriods
# and their sub-regions.
def ratio(profs):
    all_centroids = calc_all_centroids(profs)  # compute all centroids
    distances = quantify_distances(all_centroids, profs)  # compute distances of group members to centorids
    avg_dist_sub_group = sum([distances[key] for key in distances]) / len(distances)  # average distance of sub groups with centriods
    avg_centroid_dist = avg_centroid(all_centroids)  # average distance between all centroids
    return avg_centroid_dist / avg_dist_sub_group


def calc_p(std, shuffled):
    n = 0
    d = 0
    for i in shuffled:
        if i > std:
            n += 1
        else:
            d += 1
    return n / d


def main():

    # get starting data, and get sub-groups by position
    all_profiles = get_all_static_profiles()
    all_profiles = all_profiles.log_convert()

    for pos in constants.Positions:
        sub = all_profiles.filter(position=pos)

        # non-shuffled ratio
        std_rat = ratio(sub)

        shuffled_rat = []
        num_shuffle = 10000
        for i in range(0, num_shuffle):
            shuffled_prof = sub.shuffle()
            shuffled_rat.append(ratio(shuffled_prof))

        print('position: ' + str(pos))
        print(calc_p(std_rat, shuffled_rat))


if __name__ == '__main__':
    main()


