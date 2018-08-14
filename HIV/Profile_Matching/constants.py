from enum import Enum

# enumeration to represent regions
class Region(Enum):
    NA = 'NA'
    EU = 'EU'
    KOR = 'KOR'
    SA = 'SA'
    ECA = 'ECA'
    BR = 'BR'
    ALL = 'ALL'


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
    AE = 'AE'
    B = 'B'
    C = 'C'