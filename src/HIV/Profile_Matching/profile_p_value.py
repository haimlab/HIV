from file_parse import get_all_static_profiles, AllStaticProfiles
from numpy import asarray
from scipy.cluster.vq import kmeans
import constants
from constants import FilterProperties
from argparse import ArgumentParser


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


def clade_specificity(num_shuffle, all_profiles, positions):
    for pos in positions:
        sub = all_profiles.filter(pos)
        std_rat = ratio(sub, FilterProperties.CLADE)
        num_above = 0
        num_below = 0
        num_equal = 0
        for i in range(0, num_shuffle):
            shuffled_prof = sub.shuffle(FilterProperties.CLADE)
            r = ratio(shuffled_prof, FilterProperties.CLADE)
            if r > std_rat:
                num_above += 1
            elif r < std_rat:
                num_below += 1
            else:
                num_equal += 1
        print(f'position: {pos}, p-value: {num_below / num_above}')
        print(f'below: {num_below}, above: {num_above}, equal: {num_equal}')


def select_sub_group(raw, clade_region_pairs, positions):
    # parse input clade, region pairs
    pairs = []
    for p in clade_region_pairs:
        [c, r] = p.split(',')
        pairs.append((constants.Clade(c), constants.Region(r)))

    # combine clade, region, position and make a list of filters
    filters = []
    for p in positions:
        filters += [(p,) + i for i in pairs]

    # use filter to grab all desired sub profiles
    all_prof = AllStaticProfiles()
    for f in filters:
        sub = raw.filter(*f).get_all_profiles()
        for s in sub:
            if s.clade() == constants.Clade.C and s.position() == 295:
                continue
            if s.clade() == constants.Clade.AE and s.position() == 332:
                continue
            all_prof.add_profile(s)

    return all_prof


# input profiles should already be shuffled
def specificity_one_round(all_profiles, group_by_type, cur_prop_val):

    # calculate distance without (distance relative to other sub groups)
    all_centroids = calc_all_centroids(all_profiles, group_by_type)
    cur_centroid = all_centroids[cur_prop_val]
    del all_centroids[cur_prop_val]
    dist_without = sum([euc_dist(cur_centroid, all_centroids[c]) for c in all_centroids]) / len(all_centroids)

    # calculate distance within
    sub_profiles = all_profiles.filter(cur_prop_val)
    total = sum([euc_dist(cur_centroid, p.get_entire_distr()) for p in sub_profiles.get_all_profiles()])
    dist_within = total / len(sub_profiles.get_all_profiles())

    return dist_within / dist_without


def specificity(num_shuffle, all_profiles, group_by_type, properties_to_shuffle):

    group_by_properties = all_profiles.attr_list(group_by_type)

    for prop in group_by_properties:
        std_rat = specificity_one_round(all_profiles, group_by_type, prop)
        num_above = 0
        num_below = 0
        for _ in range(0, num_shuffle):
            shuffled_profiles = all_profiles
            for p in properties_to_shuffle:
                shuffled_profiles = shuffled_profiles.shuffle(p)
            shuffled_rat = specificity_one_round(shuffled_profiles, group_by_type, prop)
            if shuffled_rat > std_rat:
                num_above += 1
            elif shuffled_rat < std_rat:
                num_below += 1
        try:
            print(f'{prop}: p-value: {num_below / num_above}')
            print(f'times shuffled: {num_shuffle}, un-shuffled ratio: {std_rat}')
            print(f'above: {num_above}, below: {num_below}, dropped: {num_shuffle - num_above - num_below}')
        except ZeroDivisionError:
            print(f'{prop}: p-value: zero division')


def corrected_position_specificity(num_shuffle, all_profiles):
    specificity(num_shuffle, all_profiles, FilterProperties.POSITION, [FilterProperties.POSITION])


def main():
    parser = ArgumentParser()
    parser.add_argument('-t', dest='type', type=str, required=True)
    parser.add_argument('-n', dest='num_shuffle', type=int, required=True)
    parser.add_argument('-e', dest='epitope', type=str, required=True)
    parser.add_argument('-p', nargs='+', dest='clade_region_pairs', required=True)
    cmd_args = parser.parse_args()

    if cmd_args.epitope == '2g12':
        positions = constants.POS_2G12
    elif cmd_args.epitope == '2f5':
        positions = constants.POS_2F5
    else:
        raise Exception('invalid epitope type')

    all_prof = select_sub_group(get_all_static_profiles(), cmd_args.clade_region_pairs, positions)
    all_prof = all_prof.log_convert()

    if cmd_args.type == 'clade':
        clade_specificity(cmd_args.num_shuffle, all_prof, positions)
    elif cmd_args.type == 'position':
        corrected_position_specificity(cmd_args.num_shuffle, all_prof)
    # elif cmd_args.type == 'clade and region':
    #     clade_and_region_specificity(cmd_args.num_shuffle, all_prof, positions, cmd_args.clade_region_pairs)
    else:
        raise Exception('invalid input')

    # group_type = FilterProperties(cmd_args.group_type)
    # shuffle_type = [FilterProperties(p) for p in cmd_args.shuffle_type]
    #
    # specificity(cmd_args.num_shuffle, all_prof, group_type, shuffle_type)


if __name__ == '__main__':
    main()
