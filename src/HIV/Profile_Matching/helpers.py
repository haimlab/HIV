from constants import AminoAcid
from file_parse import logConvert


# compute profile from all given envelopes
# i -> column index of Amino Acids to be counted
def envelopes_to_profile(envs, i):
    count = {aa: 0 for aa in AminoAcid}
    for e in envs:
        count[AminoAcid(e[i])] += 1
    s = sum([count[aa] for aa in count])
    for aa in count:
        count[aa] = logConvert(count[aa] / s)
    return count


def envelopes_to_profile_no_log(envs, i):
    count = {aa: 0 for aa in AminoAcid}
    for e in envs:
        count[AminoAcid(e[i])] += 1
    s = sum([count[aa] for aa in count])
    for aa in count:
        count[aa] = count[aa] / s
    return count