from query import QueryInput, Query
import sys
from src.HIV import constants
from file_parse import get_all_dynamic_profiles
from helpers import calc_year
from weighted_fit import calcFit
import csv


def write_validation_results(query, out_file_name, y2predict):
    with open(out_file_name, 'w') as out_file:
        writer = csv.writer(out_file, lineterminator='\n')

        # writer input query
        writer.writerow(['input query'])
        writer.writerow(['position:', str(query.input.position)])
        writer.writerow(['clade:', query.input.clade.value])
        writer.writerow(['region:', query.input.region.value])
        writer.writerow(['period', query.input.year_range])
        writer.writerow(['Profile'] + [aa.value for aa in constants.AminoAcid])
        writer.writerow([''] + [query.input.profile[aa] for aa in constants.AminoAcid])
        writer.writerow([''])

        # write best mathcing information
        clade, region, _, _ = query.find_diff_region_best_match()
        writer.writerow(['best matching information'])
        writer.writerow(['region:', region.value])
        writer.writerow(['Profile'] + [aa.value for aa in constants.AminoAcid])
        all_prof = get_all_dynamic_profiles().filter(query.input.clade, query.input.region, query.input.position)
        for p in all_prof.get_all_profiles():
            calcFit(p)
        prof = {}
        for aa in constants.AminoAcid:
            cur_prof = all_prof.filter(aa)
            if len(cur_prof.get_all_profiles()) != 1:
                raise Exception('something is wrong')
            p = cur_prof.get_all_profiles()[0]
            prof[aa] = p.get_distr(calc_year(query.input.year_range))
        writer.writerow([''] + [prof[aa] for aa in constants.AminoAcid])
        writer.writerow([''])

        # write profile at year to be predicted using fit equations
        writer.writerow(['profile predicted'])
        writer.writerow(['year:', str(y2predict)])
        writer.writerow(['Profile'] + [aa.value for aa in constants.AminoAcid])
        all_prof = get_all_dynamic_profiles().filter(query.input.clade, query.input.region, query.input.position)
        for p in all_prof.get_all_profiles():
            calcFit(p)
        prof = {}
        for aa in constants.AminoAcid:
            cur_prof = all_prof.filter(aa)
            if len(cur_prof.get_all_profiles()) != 1:
                raise Exception('something is wrong')
            p = cur_prof.get_all_profiles()[0]
            prof[aa] = p.fit.calc_distr(y2predict)
        writer.writerow([''] + [prof[aa] for aa in constants.AminoAcid])
        writer.writerow([''])

        # write actual profile information
        writer.writerow(['measured future profile'])
        writer.writerow(['year:', str(y2predict)])
        writer.writerow(['region:', query.input.region.value])
        writer.writerow(['Profile'] + [aa.value for aa in constants.AminoAcid])
        all_prof = get_all_dynamic_profiles().filter(query.input.clade, query.input.region, query.input.position)
        for p in all_prof.get_all_profiles():
            calcFit(p)
        prof = {}
        for aa in constants.AminoAcid:
            cur_prof = all_prof.filter(aa)
            if len(cur_prof.get_all_profiles()) != 1:
                raise Exception('something is wrong')
            p = cur_prof.get_all_profiles()[0]
            prof[aa] = p.get_distr(y2predict)
        writer.writerow([''] + [prof[aa] for aa in constants.AminoAcid])
        writer.writerow([''])


def main():

    query_inputs = [
        QueryInput(constants.query_profile_B_EU_339_87_to_94, 339, '[1987, 1994]',
                   constants.Clade.B, constants.Region.EU),
        QueryInput(constants.query_profile_B_EU_392_87_to_94, 392, '[1987, 1994]',
                   constants.Clade.B, constants.Region.EU)
    ]

    for file_name, query_input in zip(sys.argv[1:], query_inputs):
        allProfiles = get_all_dynamic_profiles()
        for p in allProfiles.get_all_profiles():
            calcFit(p)
        query = Query(query_input, allProfiles)
        write_validation_results(query, file_name, 2007)


if __name__ == '__main__':
    main()
