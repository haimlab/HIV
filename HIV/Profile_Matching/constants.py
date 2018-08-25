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


query_profile_B_EU_295_00_to_04 = {  # B EU 295 2000 - 2004
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


query_profile_B_EU_332_00_to_04 = {  # B EU 295 2000 - 2004
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


query_profile_B_EU_339_00_to_04 = {  # B EU 295 2000 - 2004
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


query_profile_B_EU_392_00_to_04 = {  # B EU 295 2000 - 2004
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


query_profile_B_EU_448_00_to_04 = {  # B EU 295 2000 - 2004
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