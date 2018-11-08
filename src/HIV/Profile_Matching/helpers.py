"""
Author: Rentian Dong

Helper methods for computing clade and positional specificity
"""

from math import log10
from src.HIV.constants import AMINOACIDS


LOG_ZERO_DEFAULT = 0.1


# compute profile from all given envelopes
# i -> column index of Amino Acids to be counted
def envelopes_to_profile(envs, i, log=True):
    """

    :param envs:
    :param i:
    :param log:
    :return:
    """
    count = {aa: 0 for aa in AMINOACIDS}
    for e in envs:
        count[e[i]] += 1
    s = sum([count[aa] for aa in count])
    for aa in count:
        count[aa] = log_convert(count[aa] / s) if log else count[aa] / s
    return count


def euc_dist(p1, p2):
    if p1.keys() != p2.keys():
        raise Exception('dimension mismatch')
    p1 = [p1[key] for key in p1.keys()]
    p2 = [p2[key] for key in p2.keys()]
    return (sum([(a - b) ** 2 for a, b in zip(p1, p2)])) ** .5


def log_convert(val):
    val = log10(LOG_ZERO_DEFAULT) if val < LOG_ZERO_DEFAULT else log10(val)
    return val - log10(LOG_ZERO_DEFAULT)


# calculate year as median of the range
# assumes input string to look like "[year1, year2]" with year1 < year2
def calc_year(year_range):
    comma_ind = year_range.find(',')
    year1 = int(year_range[1:comma_ind])
    year2 = int(year_range[comma_ind + 2:-1])
    return (year1 + year2) / 2