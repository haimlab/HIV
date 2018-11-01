from argparse import ArgumentParser

AMINOACIDS = ('Z', 'N', 'T', 'S', 'D', 'E',
              'K', 'R', 'H', 'Y', 'Q', 'I',
              'L', 'V', 'A', 'C', 'F', 'G',
              'M', 'P', 'W')

def main():
    ifn = ''
    year_ranges = (1, 2), (3, 4)
    positions = [1, 2, 3]
    aa_counts = {p: {r: {aa: 0 for aa in AMINOACIDS} for r in year_ranges} for p in positions}

    with open(ifn, 'r') as f:
        fr = next(f)
        positions = filter(lambda x: any(list(map(lambda y: y[0] <= int(x) <= y[1], year_ranges))), fr)
        print(positions)
