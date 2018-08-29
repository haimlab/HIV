from enum import Enum


# enumeration to represent regions
class Region(Enum):
    ALL = 'ALL'
    ASiA = 'ASIA'
    BR = 'BR'
    CN = 'CN'
    ECA = 'ECA'
    EU = 'EU'
    INNP = 'IN&NP'
    KR = 'KR'
    NA = 'NA'
    SA = 'SA'
    TH = 'TH'
    BRAZIL = 'BRAZIL'


# enumeration to represent AA
class AminoAcid(Enum):
    Z = 'Z'
    N = 'N'
    T = 'T'
    S = 'S'
    D = 'D'
    E = 'E'
    K = 'K'
    R = 'R'
    H = 'H'
    Y = 'Y'
    Q = 'Q'
    I = 'I'
    L = 'L'
    V = 'V'
    A = 'A'
    C = 'C'
    F = 'F'
    G = 'G'
    M = 'M'
    P = 'P'
    W = 'W'


# enumeration to represent clades
class Clade(Enum):
    A1 = 'A1'
    AE = 'AE'
    B = 'B'
    C = 'C'


# properties that can be applied with filter
class FilterProperties(Enum):
    CLADE = 'CLADE'
    REGION = 'REGION'
    POSITION = 'POSITION'
    AMINOACID = 'AMINOACID'


query_profile_B_EU_295_00_to_04 = {
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


query_profile_B_EU_332_00_to_04 = {
    AminoAcid('Z'): 87.23,
    AminoAcid('N'): 1.82,
    AminoAcid('T'): 9.09,
    AminoAcid('S'): 0,
    AminoAcid('D'): 0,
    AminoAcid('E'): 0,
    AminoAcid('K'): 0,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 1.82,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_EU_339_00_to_04 = {
    AminoAcid('Z'): 67.23,
    AminoAcid('N'): 9.09,
    AminoAcid('T'): 0,
    AminoAcid('S'): 3.64,
    AminoAcid('D'): 5.45,
    AminoAcid('E'): 7.27,
    AminoAcid('K'): 0,
    AminoAcid('R'): 0,
    AminoAcid('H'): 5.45,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 1.82,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_EU_392_00_to_04 = {
    AminoAcid('Z'): 81.82,
    AminoAcid('N'): 7.23,
    AminoAcid('T'): 0,
    AminoAcid('S'): 3.64,
    AminoAcid('D'): 5.45,
    AminoAcid('E'): 0,
    AminoAcid('K'): 1.82,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_EU_448_00_to_04 = {
    AminoAcid('Z'): 89.09,
    AminoAcid('N'): 0,
    AminoAcid('T'): 1.82,
    AminoAcid('S'): 1.82,
    AminoAcid('D'): 1.82,
    AminoAcid('E'): 0,
    AminoAcid('K'): 3.64,
    AminoAcid('R'): 0,
    AminoAcid('H'): 1.82,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_295_00_to_04 = {
    AminoAcid('Z'): 83.59,
    AminoAcid('N'): 2.56,
    AminoAcid('T'): 7.18,
    AminoAcid('S'): 1.54,
    AminoAcid('D'): 1.03,
    AminoAcid('E'): 0.51,
    AminoAcid('K'): 1.54,
    AminoAcid('R'): 0.51,
    AminoAcid('H'): 0.51,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 1.03,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_332_00_to_04 = {
    AminoAcid('Z'): 83.59,
    AminoAcid('N'): 3.01,
    AminoAcid('T'): 10.26,
    AminoAcid('S'): 0.51,
    AminoAcid('D'): 0.51,
    AminoAcid('E'): 0.51,
    AminoAcid('K'): 0,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 1.03,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0.51,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_339_00_to_04 = {
    AminoAcid('Z'): 69.74,
    AminoAcid('N'): 3.08,
    AminoAcid('T'): 4.62,
    AminoAcid('S'): 1.03,
    AminoAcid('D'): 2.05,
    AminoAcid('E'): 7.18,
    AminoAcid('K'): 5.64,
    AminoAcid('R'): 2.05,
    AminoAcid('H'): 1.03,
    AminoAcid('Y'): 0.51,
    AminoAcid('Q'): 0.51,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0.51,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0.51,
    AminoAcid('M'): 1.03,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0.51,
}


query_profile_B_NA_392_00_to_04 = {
    AminoAcid('Z'): 82.55,
    AminoAcid('N'): 6.67,
    AminoAcid('T'): 2.56,
    AminoAcid('S'): 4.1,
    AminoAcid('D'): 2.05,
    AminoAcid('E'): 0,
    AminoAcid('K'): 1.03,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0.51,
    AminoAcid('Y'): 0.51,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_448_00_to_04 = {
    AminoAcid('Z'): 91.28,
    AminoAcid('N'): 0,
    AminoAcid('T'): 0,
    AminoAcid('S'): 2.56,
    AminoAcid('D'): 0,
    AminoAcid('E'): 1.03,
    AminoAcid('K'): 2.56,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0.51,
    AminoAcid('Q'): 0.51,
    AminoAcid('I'): 1.54,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_295_05_to_09 = {
    AminoAcid('Z'): 74.48,
    AminoAcid('N'): 5.34,
    AminoAcid('T'): 10.98,
    AminoAcid('S'): 0,
    AminoAcid('D'): 2.08,
    AminoAcid('E'): 1.78,
    AminoAcid('K'): 2.37,
    AminoAcid('R'): 0,
    AminoAcid('H'): 1.48,
    AminoAcid('Y'): 0.59,
    AminoAcid('Q'): 0.59,
    AminoAcid('I'): 0.30,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_332_05_to_09 = {
    AminoAcid('Z'): 83.68,
    AminoAcid('N'): 1.78,
    AminoAcid('T'): 8.01,
    AminoAcid('S'): 0.3,
    AminoAcid('D'): 0.3,
    AminoAcid('E'): 2.37,
    AminoAcid('K'): 0.59,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 2.67,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0.3,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_339_05_to_09 = {
    AminoAcid('Z'): 73,
    AminoAcid('N'): 1.48,
    AminoAcid('T'): 4.75,
    AminoAcid('S'): 1.48,
    AminoAcid('D'): 3.56,
    AminoAcid('E'): 4.45,
    AminoAcid('K'): 5.93,
    AminoAcid('R'): 1.19,
    AminoAcid('H'): 1.48,
    AminoAcid('Y'): 0.3,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0.3,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0.3,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 1.19,
    AminoAcid('M'): 0.59,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_392_05_to_09 = {
    AminoAcid('Z'): 84.87,
    AminoAcid('N'): 7.42,
    AminoAcid('T'): 0.89,
    AminoAcid('S'): 1.78,
    AminoAcid('D'): 2.67,
    AminoAcid('E'): 0,
    AminoAcid('K'): 1.48,
    AminoAcid('R'): 0.3,
    AminoAcid('H'): 0.3,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0.3,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}


query_profile_B_NA_448_05_to_09 = {
    AminoAcid('Z'): 92.58,
    AminoAcid('N'): 0,
    AminoAcid('T'): 0.59,
    AminoAcid('S'): 1.19,
    AminoAcid('D'): 1.78,
    AminoAcid('E'): 0.59,
    AminoAcid('K'): 1.78,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0.3,
    AminoAcid('Q'): 0.59,
    AminoAcid('I'): 0.3,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0.3,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}