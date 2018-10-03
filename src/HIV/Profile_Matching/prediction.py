from file_parse import get_all_dynamic_profiles
from weighted_fit import calcFit
from constants import AminoAcid, Clade, Region
import argparse
from os.path import join
from csv import writer

def predict(all_profiles, clade, region, position, year):
    res = {}
    actual = {}
    for aa in AminoAcid:
        try:
            p = all_profiles.filter(clade, region, position, aa).get_only_profile()
        except Exception as e:
            if str(e) == 'cannot find single profile':
                return None, None
            else:
                raise
        distr = p.fit.calc_distr(year)
        res[aa] = distr
        try:
            actual[aa] = p.get_distr(year)
        except:
            actual[aa] = -1
    return res, actual


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', nargs='+', dest='positions', type=int, required=True)
    parser.add_argument('-c', nargs='+', dest='clades', type=str, required=True)
    parser.add_argument('-y', nargs='+', dest='years', type=float, required=True)
    parser.add_argument('-r', nargs='+', dest='regions', type=str, required=True)
    parser.add_argument('-t', dest='omit_last', type=bool, required=True)
    parser.add_argument('-o', dest='out_dir', type=str, required=True)
    cmd_args = parser.parse_args()
    positions = cmd_args.positions
    clades = [Clade(c) for c in cmd_args.clades]
    years = cmd_args.years
    regions = [Region(r) for r in cmd_args.regions]
    out_dir = cmd_args.out_dir

    all_profiles = get_all_dynamic_profiles()
    for p in all_profiles.get_all_profiles():
        if cmd_args.omit_last:
            calcFit(p, skip_tail=True)
        else:
            calcFit(p, skip_tail=False)

    rows_predicted = []
    rows_actual = []
    header = ['Clade', 'Region', 'Position', 'Year'] + [str(aa.value) for aa in AminoAcid]
    for c in clades:
        for r in regions:
            for p in positions:
                for y in years:
                    row = [str(c.value), str(r.value), str(p), str(y)]
                    res, actual = predict(all_profiles, c, r, p, y)
                    if res is None and actual is None:
                        continue
                    rows_predicted.append(row + [str(res[aa]) for aa in AminoAcid])
                    rows_actual.append(row + [str(actual[aa]) for aa in AminoAcid])
    with open(join(out_dir, 'out.csv'), 'w') as out_file:
        w = writer(out_file, lineterminator='\n')
        w.writerow(header)
        w.writerow(['predicted'])
        w.writerows(rows_predicted)
        w.writerow(['measured'])
        w.writerows(rows_actual)


if __name__ == '__main__':
    main()