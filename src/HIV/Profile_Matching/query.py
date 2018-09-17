import sys
import csv
import constants
from file_parse import get_all_dynamic_profiles, calcYear
from weighted_fit import calcFit
from os.path import join
from math import log10


class QueryInput:
    def __init__(self, profile, position, year_range, clade, region):
        self.profile = profile
        self.position = position
        self.year_range = year_range
        self.clade = clade
        self.region = region


class Query:
    def __init__(self, query_input, all_profiles):
        self.all_profiles = all_profiles.filter(position=query_input.position)
        self.input = query_input
        self.results = {}  # {(clade, region) -> {amino acid -> year}}
        self.scores = {}  # {(clade, region) -> stdev of corresponding profile), for writing to file specifically
        self.best_match_years = {}  # {(clade, region) -> year)}
        self.fits = {}  # {(clade, region) -> {amino acid -> fit object}}
        self.best_fit_clade = None
        self.best_fit_region = None

        for aminoAcid in self.input.profile:
            sub_profiles = self.all_profiles.filter(aminoAcid=aminoAcid)
            for p in sub_profiles.get_all_profiles():
                if p.fit is None:
                    print(p.tag())
                self.add_result(p.clade(), p.region(), aminoAcid,
                                p.fit.calcYear(query_input.profile[aminoAcid],
                                calcYear(self.input.year_range)), p.fit)
        self.find_best_match()

    def add_result(self, clade, region, aminoAcid, year, fit):
        try:
            self.results[(clade, region)][aminoAcid] = year
            self.fits[(clade, region)][aminoAcid] = fit
        except KeyError:
            self.results[(clade, region)] = {aminoAcid: year}
            self.fits[(clade, region)] = {aminoAcid: fit}

    def __get_corrected_r2(self, clade, region, amino_acid):
        fit = self.fits[(clade, region)][amino_acid]
        query_percentage = self.input.profile[amino_acid]
        if query_percentage != 0 and fit.slope == 0:
            return .1
        else:
            return fit.r

    def get_avg_percent(self, p1, p2):
        res = (p1 + p2) / 2
        if res == 0:
            res = .5
        return res

    def calc_avg_percentage(self, aa, clade, region):
        cand_prof = self.all_profiles.filter(clade=clade, region=region, aminoAcid=aa)
        if len(cand_prof.get_all_profiles()) != 1:
            raise Exception('something went wrong')
        only_prof = cand_prof.get_all_profiles()[0]
        return sum(only_prof.distr) / len(only_prof.distr)

    def __calc_log_avgs(self, clade, region, prof):
        _min = float('inf')
        avgs = {}
        for amino_acid in prof:
            avg_percentage = self.calc_avg_percentage(amino_acid, clade, region)
            if avg_percentage > 100:
                print(avg_percentage)
            percent = self.input.profile[amino_acid]
            log_avg = log10(self.get_avg_percent(percent, avg_percentage))
            if _min > log_avg:
                _min = log_avg
            avgs[amino_acid] = log_avg
        normalized_avgs = {aa: avgs[aa] - _min for aa in avgs}
        return normalized_avgs

    # compute all weights corresponding with a clade, region, and profile
    # profile: {amino acid: year calculated through best fit based on query percentage}
    def __calc_weights(self, clade, region, prof):
        r2s = {aa: self.__get_corrected_r2(clade, region, aa) for aa in prof}
        log_percent_avgs = self.__calc_log_avgs(clade, region, prof)
        normalized_weights = {aa: r2s[aa] * log_percent_avgs[aa] for aa in prof}
        return normalized_weights

    def clac_r2_weighted_stdev(self, prof, clade, region, weighted_avg_year):
        weighted_avg = weighted_avg_year

        # then calculate weighted stdev
        weighted_square_sum = 0
        weights = self.__calc_weights(clade, region, prof)
        for aa in prof:
            if abs(prof[aa]) == float('inf'):  # skip inf
                continue
            weighted_square_sum += weights[aa] * (prof[aa] - weighted_avg) ** 2
        return (weighted_square_sum / (len(prof) - 1)) ** .5

    def find_best_match(self):

        def calc_best_matching_year(profile, clade, region):

            # weighted sum and number of year entries
            weights = self.__calc_weights(clade, region, profile)
            weighted_year_sum = 0
            weighted_year_num = 0
            for aa in profile:
                if abs(profile[aa]) == float('inf'):  # skip inf
                    continue
                weighted_year_sum += weights[aa] * profile[aa]
                weighted_year_num += weights[aa]

            return weighted_year_sum / weighted_year_num

        score = sys.float_info.max
        for clade, region in self.results:
            profileDict = self.results[(clade, region)]
            # cur_score = calc_stdev(profileDict)
            cur_best_match_year = calc_best_matching_year(profileDict, clade, region)
            cur_score = self.clac_r2_weighted_stdev(profileDict, clade, region, cur_best_match_year)
            self.scores[(clade, region)] = cur_score
            self.best_match_years[(clade, region)] = cur_best_match_year
            if cur_score < score:
                score = cur_score
                self.best_fit_clade = clade
                self.best_fit_region = region
        return self.best_fit_clade, self.best_fit_region

    # write intermediate results for verfication purposes
    def write_intermediate_results(self, fileName):
        all_rows = []
        amino_acids_list = []
        for aminoAcid in constants.AminoAcid:
            amino_acids_list.append(aminoAcid.value)

        # query rows
        all_rows.append(['position:', str(self.input.position)])
        all_rows.append(['clade:', self.input.clade.value])
        all_rows.append(['region:', self.input.region.value])
        all_rows.append(['period', self.input.year_range])
        all_rows.append(['', '', '', ''] + amino_acids_list)
        query_row = ['query', '', '', '']
        for aminoAcid in constants.AminoAcid:
            query_row.append(self.input.profile[aminoAcid])
        all_rows.append(query_row)
        all_rows.append([])

        # calculated year, stdev, and fit parameters
        all_rows.append(['', 'stdev', 'best match year', ''] + amino_acids_list)
        print(len(self.results))
        for clade, region in self.results:
            row_collection = [
                [clade.value + ", " + region.value, self.scores[clade, region], self.best_match_years[clade, region], ''],
                ['avg_percentage', '', '', ''],
                ['slope', '', '', ''],
                ['y_intercept', '', '', ''],
                ['r', '', '', ''],
                []
            ]
            avgs = self.__calc_log_avgs(clade, region, [a for a in constants.AminoAcid])
            for aminoAcid in constants.AminoAcid:
                row_collection[0].append(self.results[(clade, region)][aminoAcid])
                row_collection[1].append(avgs[aminoAcid])
                row_collection[2].append(self.fits[(clade, region)][aminoAcid].slope)
                row_collection[3].append(self.fits[(clade, region)][aminoAcid].y_intercept)
                row_collection[4].append(self.fits[(clade, region)][aminoAcid].r)
            for row in row_collection:
                all_rows.append(row)
            all_rows.append([])

        # write all other profiles that shared the same period with query
        period = calcYear(self.input.year_range)
        all_other_profiles = {}  # {(clade, region) -> {AminoAcid -> percentage}}
        all_profiles = get_all_dynamic_profiles()
        all_profiles = all_profiles.filter(position=self.input.position)
        for profile in all_profiles.get_all_profiles():
            key = profile.clade(), profile.region()
            amino_acid = profile.amino_acid()
            percentage = None
            for i in range(0, len(profile.years)):
                if profile.years[i] == period:
                    percentage = profile.distr[i]
                    break
            try:
                all_other_profiles[key][amino_acid] = percentage
            except KeyError:
                all_other_profiles[key] = {amino_acid: percentage}

        # put all other profiles into the list of rows to be writteqn
        all_rows.append([])
        all_rows.append(['other profiles in the same period'])
        all_rows.append(['', '', ''] + amino_acids_list)
        for clade, region in all_other_profiles:
            profile = all_other_profiles[(clade, region)]
            percentage_row = [clade.value + ", " + region.value, '']
            for amino_acid in constants.AminoAcid:
                percentage_row.append(profile[amino_acid])
            all_rows.append(percentage_row)

        # write to output file
        with open(fileName, 'w') as outFile:
            writer = csv.writer(outFile, lineterminator='\n')
            writer.writerows(all_rows)


if __name__ == '__main__':

    query_inputs = [
        QueryInput(constants.query_profile_B_NA_295_05_to_09, 295, '[2005, 2009]', constants.Clade.B, constants.Region.NA),
        QueryInput(constants.query_profile_B_NA_332_05_to_09, 332, '[2005, 2009]', constants.Clade.B, constants.Region.NA),
        QueryInput(constants.query_profile_B_NA_339_05_to_09, 339, '[2005, 2009]', constants.Clade.B, constants.Region.NA),
        QueryInput(constants.query_profile_B_NA_392_05_to_09, 392, '[2005, 2009]', constants.Clade.B, constants.Region.NA),
        QueryInput(constants.query_profile_B_NA_448_05_to_09, 448, '[2005, 2009]', constants.Clade.B, constants.Region.NA)
    ]

    for query_input in query_inputs:
        allProfiles = get_all_dynamic_profiles()
        for p in allProfiles.get_all_profiles():
            calcFit(p)
        query = Query(query_input, allProfiles)
        print(query.best_fit_clade)
        print(query.best_fit_region)
        file_name = join(sys.argv[1], str(query_input.position) + '_r2=0.4.csv')
        query.write_intermediate_results(file_name)
