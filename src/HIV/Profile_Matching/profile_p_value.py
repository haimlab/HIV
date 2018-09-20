from file_parse import get_all_static_profiles
from numpy import asarray
from scipy.cluster.vq import kmeans
import constants
from constants import FilterProperties


# output type {prop_type - {amino acid -> percentage}
def calc_all_centroids(all_profiles, prop_type):
    centroids = {}
    props = all_profiles.attr_list(prop_type)
    for p in props:
        sub = all_profiles.filter(p)
        centroid = calc_centroid(sub)
        centroids[p] = centroid
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
    for prop in centroids:
        sub_profiles = all_profiles.filter(prop)
        centroid = centroids[prop]
        for prof in sub_profiles.get_all_profiles():
            distr = prof.get_entire_distr()
            distances[(prop, prof.region())] = euc_dist(centroid, distr)
    return distances


# average distance between all centroids
# input {clade -> {aa -> centroid_val}}
def avg_centroid(centroids):
    c_list = [centroids[prop] for prop in centroids]
    dist_list = []
    for i in range(0, len(c_list)):
        c = c_list[i]
        for j in range(0, len(c_list)):
            if j != i:
                dist_list.append(euc_dist(c, c_list[j]))
    return sum(dist_list) / len(dist_list)


# calculates the ratio of average distance between all centroids and average distance between all centriods
# and their sub-regions.
def ratio(profs, prop_type):
    all_centroids = calc_all_centroids(profs, prop_type)  # compute all centroids
    distances = quantify_distances(all_centroids, profs)  # compute distances of group members to centorids
    avg_dist_sub_group = sum([distances[key] for key in distances]) / len(distances)
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


def clade_specificity(num_shuffle):

    # get starting data, and get sub-groups by position
    all_profiles = get_all_static_profiles()
    all_profiles = all_profiles.log_convert()

    for pos in constants.Positions:
        sub = all_profiles.filter(position=pos)

        # non-shuffled ratio
        std_rat = ratio(sub, FilterProperties.CLADE)

        shuffled_rat = []
        for i in range(0, num_shuffle):
            shuffled_prof = sub.shuffle(FilterProperties.CLADE)
            shuffled_rat.append(ratio(shuffled_prof, FilterProperties.CLADE))

        print('position: ' + str(pos))
        print(calc_p(std_rat, shuffled_rat))


def pos_specificity(num_shuffle):

    # get starting data, and get sub-groups by position
    all_profiles = get_all_static_profiles()
    all_profiles = all_profiles.log_convert()

    # non-shuffled ratio
    std_rat = ratio(all_profiles, FilterProperties.POSITION)

    shuffled_rat = []
    for i in range(0, num_shuffle):
        shuffled_prof = all_profiles.shuffle(FilterProperties.POSITION)
        shuffled_rat.append(ratio(shuffled_prof, FilterProperties.POSITION))

    print('p-value: ' + str(calc_p(std_rat, shuffled_rat)))


if __name__ == '__main__':
    clade_specificity(10000)
