import sys
from constants import AminoAcid
from weighted_fit import get_all_profiles


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


class QueryResult:
    def __init__(self, position):
        self.position = position
        self.results = {}  # {(clade, country) -> {amino acid -> year}}

    def add_result(self, clade, region, aminoAcid, year):
        try:
            self.results[(clade, region)][aminoAcid] = year
        except KeyError:
            self.results[(clade, region)] = {aminoAcid: year}

    def calc_stdev_(self, profile):

        # calculate average
        sum = 0
        for aminoAcid in profile:
            sum += profile[aminoAcid]
        mean = sum / len(profile)

        # calculate standard deviation
        sum = 0
        for aminoAcid in profile:
            sum += (profile[aminoAcid] - mean) ** 2
        return (sum / len(profile)) ** .5

    def find_best_match(self):
        best_fit_clade = None
        best_fit_region = None
        score = sys.float_info.max
        for clade, region in self.results:
            profileDict = self.results[(clade, region)]
            cur_score = self.calc_stdev_(profileDict)
            if cur_score < score:
                score = cur_score
                best_fit_clade = clade
                best_fit_region = region
        return best_fit_clade, best_fit_region


def query(queryProfile, allProfiles):
    result = QueryResult(queryProfile.position)
    allProfiles = allProfiles.filter(position=queryProfile.position)
    for aminoAcid in queryProfile.profile:
        allProfiles = allProfiles.filter(aminoAcid=aminoAcid)
        for p in allProfiles.profiles:
            result.add_result(p.clade, p.region, aminoAcid, p.fit.calcYear(queryProfile.profile[aminoAcid]))
    return result.find_best_match()


if __name__ == '__main__':
    allProfiles = get_all_profiles('C:\\Users\\rdong6\\Desktop\\Weighted Fit Test Data\\Processed')
    result = query(queryInput, allProfiles)
    print(result)