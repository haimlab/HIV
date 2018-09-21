from query import QueryInput, Query
import sys
import constants
from file_parse import get_all_dynamic_profiles, calcYear
from weighted_fit import calcFit
from os.path import join
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
        clade, region = query.find_diff_region_best_match()
        writer.writerow(['best matching information'])
        writer.writerow(['region:', region.value])
        writer.writerow(['Profile'] + [aa.value for aa in constants.AminoAcid])
        all_prof = get_all_dynamic_profiles().filter(clade=query.input.clade, region=query.input.region,
                                                     position=query.input.position).get_all_profiles()
        prof = {}
        for aa in constants.AminoAcid:
            cur_prof = all_prof.filter(aminoAcid=aa)
            if len(cur_prof.get_all_profiles()):
                raise Exception('something is wrong')
            p = cur_prof.get_all_profiles()[0]
            prof[aa] = p.get_distr(calcYear(query.input.year_range))
        writer.writerow([''] + [prof[aa] for aa in constants.AminoAcid])
        writer.writerow([''])

        # write profile at year to be predicted using fit equations
        writer.writerow(['profile predicted'])
        writer.writerow(['year:', str(y2predict)])
        writer.writerow(['Profile'] + [aa.value for aa in constants.AminoAcid])
        all_prof = get_all_dynamic_profiles().filter(clade=query.input.clade, region=query.input.region,
                                                     position=query.input.position).get_all_profiles()
        prof = {}
        for aa in constants.AminoAcid:
            cur_prof = all_prof.filter(aminoAcid=aa)
            if len(cur_prof.get_all_profiles()):
                raise Exception('something is wrong')
            p = cur_prof.get_all_profiles()[0]
            prof[aa] = p.fit.calc_distr(y2predict)
        writer.writerow([''] + [prof[aa] for aa in constants.AminoAcid])
        writer.writerow([''])

        # write actual profile information
        # TODO continue here


def main():

    query_inputs = [
        QueryInput(constants.query_profile_B_EU_339_87_to_94, 339, '[1987, 1994]', constants.Clade.B, constants.Region.EU),
        QueryInput(constants.query_profile_B_EU_392_87_to_94, 392, '[1987, 1994]', constants.Clade.B, constants.Region.EU)
    ]

    for query_input in query_inputs:
        allProfiles = get_all_dynamic_profiles()
        for p in allProfiles.get_all_profiles():
            calcFit(p)
        query = Query(query_input, allProfiles)
        clade, region = query.find_diff_region_best_match()
        print(clade)
        print(region)
        file_name = join(sys.argv[1], str(query_input.position) + '_r2=0.4.csv')
        query.write_intermediate_results(file_name)

if __name__ == '__main__':
    pass