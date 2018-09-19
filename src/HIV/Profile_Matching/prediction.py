from file_parse import get_all_dynamic_profiles
from weighted_fit import calcFit
from constants import AminoAcid, Clade, Region


def predict(all_profiles, clade, region, position, year):
    res = {}
    for aa in AminoAcid:
        p = all_profiles.filter(clade=clade, region=region, position=position, aminoAcid=aa)
        if len(p.get_all_profiles()) != 1:
            raise Exception('something is wrong')
        distr = p.get_all_profiles()[0].fit.calc_distr(year)
        res[aa] = distr
    return res


def main():
    all_profiles = get_all_dynamic_profiles()
    for p in all_profiles.get_all_profiles():
        calcFit(p)

    years = [2018, 2028, 2038, 2048]
    clades = [Clade.AE, Clade.B, Clade.C]
    positions = [295, 332, 339, 392, 448]
    for c in clades:
        for pos in positions:
            for y in years:
                res = predict(all_profiles, c, Region.ALL, pos, y)
                print(c)
                print(pos)
                print(y)
                print(res)
        print()
        print()


if __name__ == '__main__':
    main()