import sys
import csv
from constants import AminoAcid
from weighted_fit import get_all_profiles


# script level constants
HEADER_ROW = ['', 'stdev', 'Z', 'N', 'T', 'S', 'D', 'E', 'K', 'R', 'H', 'Y',
                  'Q', 'I', 'L', 'V', 'A', 'C', 'F', 'G', 'M', 'P', 'W']


class QueryInput:
    def __init__(self, position, profile):
        self.position = position
        self.profile = profile  # mapping from Amino Acid to percentage


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
queryInput = QueryInput(295, queryProfile)

# TODO, refactor QueryResult and QueryInput into one class
class QueryResult:
    def __init__(self, position, input):
        self.input = input
        self.position = position
        self.results = {}  # {(clade, region) -> {amino acid -> (year, fit object)}}
        self.scores = {}  # {(clade, region) -> stdev of corresponding profile), for writing to file specifically
        self.best_fit_clade = None
        self.best_fit_region = None

    def add_result(self, clade, region, aminoAcid, year, fit):
        try:
            self.results[(clade, region)][aminoAcid] = year, fit
        except KeyError:
            self.results[(clade, region)] = {aminoAcid: (year, fit)}

    def calc_stdev(self, profile):

        # calculate average
        sum = 0
        for aminoAcid in profile:
            sum += profile[aminoAcid][0]
        mean = sum / len(profile)

        # calculate standard deviation
        sum = 0
        for aminoAcid in profile:
            sum += (profile[aminoAcid][0] - mean) ** 2
        return (sum / len(profile)) ** .5

    def find_best_match(self):
        score = sys.float_info.max
        for clade, region in self.results:
            profileDict = self.results[(clade, region)]
            cur_score = self.calc_stdev(profileDict)
            self.scores[(clade, region)] = cur_score
            if cur_score < score:
                score = cur_score
                self.best_fit_clade = clade
                self.best_fit_region = region
        return self.best_fit_clade, self.best_fit_region


    # write intermediate results for verfication purposes
    def write_intermediate_results(self, fileName):
        with open(fileName, 'w') as outFile:
            writer = csv.writer(outFile, lineterminator='\n')

            # write header
            writer.writerow(HEADER_ROW)

            # write query
            queryRow = ['query', '']
            for aminoAcid in AminoAcid:
                queryRow.append(self.input.profile[aminoAcid])
            writer.writerow(queryRow)

            # writer calculated year and stdev using fit equations
            for clade, region in self.results:

                # write years calculated using fit
                years = [clade.value + ", " + region.value, self.scores[clade, region]]
                slopes = ['slope', '']
                y_intercepts = ['y_intercept', '']
                r_vals = ['r', '']
                profile = self.results[(clade, region)]
                for aminoAcid in AminoAcid:
                    years.append(profile[aminoAcid][0])
                    slopes.append(profile[aminoAcid][1].slope)
                    y_intercepts.append(profile[aminoAcid][1].y_intercept)
                    r_vals.append(profile[aminoAcid][1].r)
                writer.writerow(years)
                writer.writerow(slopes)
                writer.writerow(y_intercepts)
                writer.writerow(r_vals)
                writer.writerow([])

def query(queryProfile, allProfiles):
    result = QueryResult(queryProfile.position, queryProfile)
    allProfiles = allProfiles.filter(position=queryProfile.position)
    for aminoAcid in queryProfile.profile:
        subProfiles = allProfiles.filter(aminoAcid=aminoAcid)
        for p in subProfiles.profiles:
            result.add_result(p.clade, p.region, aminoAcid, p.fit.calcYear(queryProfile.profile[aminoAcid]), p.fit)
    result.find_best_match()
    return result


if __name__ == '__main__':
    allProfiles = get_all_profiles('C:\\Users\\rdong6\\Desktop\\Weighted Fit Test Data\\Processed')
    result = query(queryInput, allProfiles)
    print(result.best_fit_clade)
    print(result.best_fit_region)
    result.write_intermediate_results('C:\\Users\\rdong6\\Desktop\\intermediate_results.csv')