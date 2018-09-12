################################################
# Program to compute the most common sequence of
# amino acids among a set of strains from a given
# location/season. Takes a fasta as input.
################################################
#Author Kallin Khan, Department of Microbiology
#Updated: 4/3/17
################################################

import collections

def ComputeMostCommon(fasta, Season):
    newFasta = []
    colEls = ['']*566
    for row in fasta:
        newFasta.append(row)
        if row[0] == '>': ##indicates the start of a new strand
            if row[-6:-1] == Season or row[-4:-1] == (Season[:2] + '|'): ##checks that the strand is in the correct season
                colNum = 1
                goodSeason = 0
            else:
                goodSeason = 1
                print(row)
        elif goodSeason == 0: ##if the strand is in the correct season, this loop will execute                
            for el in row:
                if el != "-" and el != '\n':
                    colNum += 1
                    colEls[colNum-2] += el
        
    correctSequence = ""
    for el in colEls:
        correctSequence += ((collections.Counter(el).most_common(1)[0])[0])



    return (correctSequence, newFasta)
