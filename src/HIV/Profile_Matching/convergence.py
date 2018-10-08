from numpy import array
from numpy.linalg import norm
from argparse import ArgumentParser
from file_parse import get_all_dynamic_profiles
from constants import AminoAcid
from csv import writer

DEST_YR = '[2010 - 2015]'
YR_RANGE = [
    '[1979 - 1986]',
    '[1987 - 1994]',
    '[1995 - 1999]',
    '[2000 - 2004]',
    '[2005 - 2009]',
    '[2010 - 2015]'
]


def calc_dist(v1, v2):
    return norm(array(v1) - array(v2))


def main():
    parser = ArgumentParser()
    parser.add_argument('-c', dest='clade', type=str, required=True)
    parser.add_argument('-p', dest='positions', type=int, required=True)
    parser.add_argument('-src_r', nargs='+', dest='src_regions', required=True)
    parser.add_argument('-dst_r', dest='dst_region', required=True)
    parser.add_argument('-o', dest='out_file', required=True)
    cmd_args = parser.parse_args()
    all = get_all_dynamic_profiles()

    res = {}
    for pos in cmd_args.positions:
        dst_prof = all.get_profile(cmd_args.clade, cmd_args.dst_region, pos, DEST_YR)
        for r in cmd_args.dst_region:
            for y in YR_RANGE:
                src_prof = all.get_profile(cmd_args.clade, r, pos, y)
                res[(cmd_args.clade, r, pos, y)] = calc_dist(src_prof, dst_prof)

    with open(cmd_args.out_file, 'w') as out_file:
        w = writer(out_file, lineterminator='\n')
        w.writerows([list(k) + res[k] for k in res])
