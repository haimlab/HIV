from file_parse import get_all_dynamic_profiles
from weighted_fit import calcFit
from constants import AminoAcid, Clade, Region
import argparse
from os.path import join
from csv import writer

def predict(all_profiles, clade, region, position, year):
    res = {}
    for aa in AminoAcid:
        p = all_profiles.filter(clade=clade, region=region, position=position, aminoAcid=aa).get_only_profile()
        distr = p.fit.calc_distr(year)
        res[aa] = distr
    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', nargs='+', dest='positions', type=int, required=True)
    parser.add_argument('-c', nargs='+', dest='clades', type=str, required=True)
    parser.add_argument('-y', nargs='+', dest='years', type=int, required=True)
    parser.add_argument('-r', nargs='+', dest='regions', type=str, required=True)
    parser.add_argument('-o', dest='out_dir', type=str, required=True)
    cmd_args = parser.parse_args()
    positions = cmd_args.positions
    clades = [Clade(c) for c in cmd_args.clades]
    years = cmd_args.years
    regions = [Region(r) for r in cmd_args.regions]
    out_dir = cmd_args.out_dir

    all_profiles = get_all_dynamic_profiles()
    for p in all_profiles.get_all_profiles():
        calcFit(p)

    rows = []
    rows.append(['Clade', 'Region', 'Position', 'Year'] + [str(aa.value) for aa in AminoAcid])
    for c in clades:
        for r in regions:
            for p in positions:
                for y in years:
                    row = [str(c.value), str(r.value), str(p), str(y)]
                    res = predict(all_profiles, c, r, p, y)
                    row += [str(res[aa]) for aa in AminoAcid]
                    rows.append(row)

    with open(join(out_dir, 'out.csv'), 'w') as out_file:
        w = writer(out_file, lineterminator='\n')
        w.writerows(rows)


if __name__ == '__main__':
    main()