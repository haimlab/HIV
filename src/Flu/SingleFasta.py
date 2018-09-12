################################################
# Program to read/edit FASTA files for individual
# seasons and draw branches based on nucleotide
# sequences and create protein branched excel sheets
# both with and without duplicates
################################################
#Author Kallin Khan, Department of Microbiology
#Updated: 4/3/17
################################################

from openpyxl import Workbook
from ComputeBranches import ComputeBranches
from DeleteDuplicates import DeleteDuplicates
from ComputeMostCommon import ComputeMostCommon
from ComputeStandardDevs import ComputeStandardDevs, aveSDs, getVal
from AnalyzeSequences import AnalyzeNucleotides, CreateExcelSheetsForProteins, saveExcel, findZSites

import csv

    
    
book = Workbook()

Season = input("Season (Ex: 14_15): ")
State = input("State (Ex: CA,FL,MA,NY,TX,WI): ")

Season1 = Season[:2] + '-' + Season[3:]
years = '20' + Season[:2] + '-20' + Season[3:]


fastaPro = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Influenza\\H3N2 Analyses\\USA\\" + State+ "\\"+ years +'\\' + State + Season1 + '-PRO.fasta', "r")
fastaNuc = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Influenza\\H3N2 Analyses\\USA\\" + State+ "\\"+ years +'\\' + State + Season1 + '-NUC.fasta', "r")
Reader = csv.reader(open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Influenza\\H3N2 Analyses\\USA\\" + State+ "\\"+ years +'\\Distance Data-Nuc', newline=''), delimiter=' ', quotechar='|')
file = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Influenza\\H3N2 Analyses\\USA\\" + State+ "\\"+ years +'\\' + State + Season1 + '-Pro-NoDup.fasta', 'w')


(rowNum, Nucs, Names) = AnalyzeNucleotides(Season, fastaNuc)
(rowNumbers, numBranches) = ComputeBranches(Reader, rowNum, 2)
(correctSeq, fasta2) = ComputeMostCommon(fastaPro, Season)
zSites = findZSites(correctSeq)
(rowNum, Xs) = CreateExcelSheetsForProteins(book,Season, correctSeq, fasta2, rowNumbers, zSites, 2)
end = ComputeStandardDevs(book, 2, 2, years, State)
DeleteDuplicates(book, rowNum, correctSeq, numBranches,file)
aveSDs(book, end)
saveExcel(book, Season, State,years)

for x in Xs:
    print(x)

