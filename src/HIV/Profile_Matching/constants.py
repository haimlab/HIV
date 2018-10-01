from enum import Enum

# 2G12 and 2F5 positions
Positions = [295, 332, 339, 392, 448]


# enumeration to represent regions
class Region(Enum):
    ALL = 'ALL'
    ASIA = 'ASIA'
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
    ALL_WO_KR_NA = 'ALLWOKRNA'


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

query_profile_C_SA_295_05_to_09 = {
    AminoAcid('Z'): 23.74,
    AminoAcid('N'): 0.46,
    AminoAcid('T'): 5.93,
    AminoAcid('S'): 0,
    AminoAcid('D'): 0,
    AminoAcid('E'): 9.13,
    AminoAcid('K'): 1.37,
    AminoAcid('R'): 0.45,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0.457,
    AminoAcid('I'): 1.83,
    AminoAcid('L'): 0.91,
    AminoAcid('V'): 52.05,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 3.653,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}

query_profile_C_SA_332_05_to_09 = {
    AminoAcid('Z'): 75.799,
    AminoAcid('N'): 10.05,
    AminoAcid('T'): 4.11,
    AminoAcid('S'): 2.28,
    AminoAcid('D'): 1.37,
    AminoAcid('E'): 1.37,
    AminoAcid('K'): 2.28,
    AminoAcid('R'): 0.457,
    AminoAcid('H'): 0.457,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0.457,
    AminoAcid('I'): 1.37,
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

query_profile_C_SA_339_05_to_09 = {
    AminoAcid('Z'): 68.037,
    AminoAcid('N'): 3.653,
    AminoAcid('T'): 11.42,
    AminoAcid('S'): 4.11,
    AminoAcid('D'): 3.196,
    AminoAcid('E'): 3.196,
    AminoAcid('K'): 0.913,
    AminoAcid('R'): 0.913,
    AminoAcid('H'): 0.457,
    AminoAcid('Y'): 1.37,
    AminoAcid('Q'): 0.457,
    AminoAcid('I'): 0.457,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0.457,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 1.37,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}

query_profile_C_SA_392_05_to_09 = {
    AminoAcid('Z'): 68.037,
    AminoAcid('N'): 14.612,
    AminoAcid('T'): 4.11,
    AminoAcid('S'): 2.74,
    AminoAcid('D'): 0.913,
    AminoAcid('E'): 0.913,
    AminoAcid('K'): 3.653,
    AminoAcid('R'): 0,
    AminoAcid('H'): 1.826,
    AminoAcid('Y'): 0.457,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0,
    AminoAcid('L'): 0.913,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0.457,
    AminoAcid('F'): 0.457,
    AminoAcid('G'): 0.457,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0.457,
    AminoAcid('W'): 0,
}

query_profile_C_SA_448_05_to_09 = {
    AminoAcid('Z'): 70.776,
    AminoAcid('N'): 0,
    AminoAcid('T'): 3.653,
    AminoAcid('S'): 19.63,
    AminoAcid('D'): 1.826,
    AminoAcid('E'): 0.913,
    AminoAcid('K'): 1.826,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0.913,
    AminoAcid('L'): 0,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0.457,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}

query_profile_C_SA_295_00_to_04 = {
    AminoAcid('Z'): 13.55,
    AminoAcid('N'): 0.733,
    AminoAcid('T'): 5.86,
    AminoAcid('S'): 0,
    AminoAcid('D'): 0,
    AminoAcid('E'): 9.524,
    AminoAcid('K'): 2.198,
    AminoAcid('R'): 0.733,
    AminoAcid('H'): 0.366,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 1.831,
    AminoAcid('L'): 0.733,
    AminoAcid('V'): 59.707,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 4.762,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}

query_profile_C_SA_332_00_to_04 = {
    AminoAcid('Z'): 78.388,
    AminoAcid('N'): 10.256,
    AminoAcid('T'): 2.564,
    AminoAcid('S'): 3.66,
    AminoAcid('D'): 1.099,
    AminoAcid('E'): 1.099,
    AminoAcid('K'): 1.099,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0.366,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0.733,
    AminoAcid('L'): 0.733,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}

query_profile_C_SA_339_00_to_04 = {
    AminoAcid('Z'): 75.092,
    AminoAcid('N'): 5.128,
    AminoAcid('T'): 7.326,
    AminoAcid('S'): 2.93,
    AminoAcid('D'): 4.396,
    AminoAcid('E'): 1.099,
    AminoAcid('K'): 1.465,
    AminoAcid('R'): 1.099,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0.366,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 1.099,
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

query_profile_C_SA_392_00_to_04 = {
    AminoAcid('Z'): 68.498,
    AminoAcid('N'): 14.286,
    AminoAcid('T'): 3.663,
    AminoAcid('S'): 7.326,
    AminoAcid('D'): 1.465,
    AminoAcid('E'): 0.366,
    AminoAcid('K'): 0.733,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0.366,
    AminoAcid('Y'): 0.733,
    AminoAcid('Q'): 1.099,
    AminoAcid('I'): 0.366,
    AminoAcid('L'): 0.366,
    AminoAcid('V'): 0,
    AminoAcid('A'): 0.733,
    AminoAcid('C'): 0,
    AminoAcid('F'): 0,
    AminoAcid('G'): 0,
    AminoAcid('M'): 0,
    AminoAcid('P'): 0,
    AminoAcid('W'): 0,
}

query_profile_C_SA_448_00_to_04 = {
    AminoAcid('Z'): 67.766,
    AminoAcid('N'): 0.366,
    AminoAcid('T'): 3.663,
    AminoAcid('S'): 24.542,
    AminoAcid('D'): 0,
    AminoAcid('E'): 0.733,
    AminoAcid('K'): 2.564,
    AminoAcid('R'): 0,
    AminoAcid('H'): 0,
    AminoAcid('Y'): 0,
    AminoAcid('Q'): 0,
    AminoAcid('I'): 0.366,
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

query_profile_B_EU_339_87_to_94 = {
    AminoAcid('Z'): 70.588,
    AminoAcid('N'): 2.9411,
    AminoAcid('T'): 8.824,
    AminoAcid('S'): 2.9411,
    AminoAcid('D'): 2.9411,
    AminoAcid('E'): 5.8824,
    AminoAcid('K'): 0,
    AminoAcid('R'): 0,
    AminoAcid('H'): 5.8824,
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

query_profile_B_EU_392_87_to_94 = {
    AminoAcid('Z'): 91.176,
    AminoAcid('N'): 5.8824,
    AminoAcid('T'): 0,
    AminoAcid('S'): 0,
    AminoAcid('D'): 2.9411,
    AminoAcid('E'): 0,
    AminoAcid('K'): 0,
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
