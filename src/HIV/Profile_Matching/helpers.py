from math import log10
from constants import AMINOACIDS
LOG_ZERO_DEFAULT = 0.1

# compute profile from all given envelopes
# i -> column index of Amino Acids to be counted
def envelopes_to_profile(envs, i, log=True):
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