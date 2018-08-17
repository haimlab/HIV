import sys
import csv
from constants import AminoAcid
from constants import Clade
from constants import Region
from weighted_fit import get_all_profiles
from weighted_fit import calcYear


# temporary hard-coded inputs
queryProfile = {
    AminoAcid('Z'): 78.5,
    AminoAcid('N'): 3.57,
    AminoAcid('T'): 5.35,
    AminoAcid('S'): 3.57,
    AminoAcid('D'): 1.78,
    AminoAcid('E'): 1.78,
    AminoAcid('K'): 3.57,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0,
    AminoAcid('V'): 1.78,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


class QueryInput:
    def __init__(self, profile, position, year_range, clade, region):
        self.profile = profile
        self.position = position
        self.year_range = year_range
        self.clade = clade
        self.region = region


class Query:
    def __init__(self, query_input, all_profiles):
        self.all_profiles = all_profiles
        self.input = query_input
        self.results = {}  # {(clade, region) -> {amino acid -> year}}
        self.scores = {}  # {(clade, region) -> stdev of corresponding profile), for writing to file specifically
        self.fits = {}  # {(clade, region) -> {amino acid -> fit object}}
        self.best_fit_clade = None
        self.best_fit_region = None

        all_profiles = all_profiles.filter(position=self.input.position)
        for aminoAcid in self.input.profile:
            sub_profiles = all_profiles.filter(aminoAcid=aminoAcid)
            for p in sub_profiles.profiles:
                self.add_result(p.clade, p.region, aminoAcid,
                                p.fit.calcYear(queryProfile[aminoAcid],
                                calcYear(self.input.year_range)), p.fit)
        self.find_best_match()

    def add_result(self, clade, region, aminoAcid, year, fit):
        try:
            self.results[(clade, region)][aminoAcid] = year
            self.fits[(clade, region)][aminoAcid] = fit
        except KeyError:
            self.results[(clade, region)] = {aminoAcid: year}
            self.fits[(clade, region)] = {aminoAcid: fit}

    def find_best_match(self):

        def calc_stdev(profile):
            sum = 0
            for aminoAcid in profile:
                sum += profile[aminoAcid]
            mean = sum / len(profile)
            sum = 0
            for aminoAcid in profile:
                sum += (profile[aminoAcid] - mean) ** 2
            return (sum / len(profile)) ** .5

        score = sys.float_info.max
        for clade, region in self.results:
            profileDict = self.results[(clade, region)]
            cur_score = calc_stdev(profileDict)
            self.scores[(clade, region)] = cur_score
            if cur_score < score:
                score = cur_score
                self.best_fit_clade = clade
                self.best_fit_region = region
        return self.best_fit_clade, self.best_fit_region


    # write intermediate results for verfication purposes
    def write_intermediate_results(self, fileName):
        all_rows = []
        amino_acids_list = []
        for aminoAcid in AminoAcid:
            amino_acids_list.append(aminoAcid.value)

        # query rows
        all_rows.append(['position:', str(self.input.position)])
        all_rows.append(['clade:', self.input.clade.value])
        all_rows.append(['region:', self.input.region.value])
        all_rows.append(['period', self.input.year_range])
        all_rows.append(['', ''] + amino_acids_list)
        query_row = ['query', '']
        for aminoAcid in AminoAcid:
            query_row.append(self.input.profile[aminoAcid])
        all_rows.append(query_row)
        all_rows.append([])

        # calculated year, stdev, and fit parameters
        all_rows.append(['', 'stdev'] + amino_acids_list)
        for clade, region in self.results:
            row_collection = [
                [clade.value + ", " + region.value, self.scores[clade, region]],
                ['slope', ''],
                ['y_intercept', ''],
                ['r', ''],
                []
            ]
            for aminoAcid in AminoAcid:
                row_collection[0].append(self.results[(clade, region)][aminoAcid])
                row_collection[1].append(self.fits[(clade, region)][aminoAcid].slope)
                row_collection[2].append(self.fits[(clade, region)][aminoAcid].y_intercept)
                row_collection[3].append(self.fits[(clade, region)][aminoAcid].r)
            for row in row_collection:
                all_rows.append(row)
            all_rows.append([])

        # write all other profiles that shared the same period with query
        period = calcYear(self.input.year_range)
        all_other_profiles = {}  # {(clade, region) -> {AminoAcid -> percentage}}
        for profile in self.all_profiles.profiles:
            key = profile.clade, profile.region
            amino_acid = profile.aminoAcid
            percentage = None
            for i in range(0, len(profile.years)):
                if profile.years[i] == period:
                    percentage = profile.distr[i]
                    break
            try:
                all_other_profiles[key][amino_acid] = percentage
            except KeyError:
                all_other_profiles[key] = {amino_acid: percentage}

        # put all other profiles into the list of rows to be written
        all_rows.append([])
        all_rows.append(['other profiles in the same period'])
        all_rows.append(['', ''] + amino_acids_list)
        for clade, region in all_other_profiles:
            profile = all_other_profiles[(clade, region)]
            percentage_row = [clade.value + ", " + region.value, '']
            for amino_acid in AminoAcid:
                percentage_row.append(profile[amino_acid])
            all_rows.append(percentage_row)

        # write to output file
        with open(fileName, 'w') as outFile:
            writer = csv.writer(outFile, lineterminator='\n')
            writer.writerows(all_rows)


if __name__ == '__main__':
    allProfiles = get_all_profiles(sys.argv[1])
    input = QueryInput(queryProfile, 295, '[2000, 2004]', Clade.B, Region.EU)
    query = Query(input, allProfiles)
    print(query.best_fit_clade)
    print(query.best_fit_region)
    query.write_intermediate_results(sys.argv[2])