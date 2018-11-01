from argparse import ArgumentParser
from csv import reader, writer


AMINOACIDS = ('Z', 'N', 'T', 'S', 'D', 'E',
              'K', 'R', 'H', 'Y', 'Q', 'I',
              'L', 'V', 'A', 'C', 'F', 'G',
              'M', 'P', 'W')
POS_START_IND = 4  # amino acid data starts at column 4
YEAR_IND = 1  # index of column holding year


def main():

    # set up cmdln arguments
    parser = ArgumentParser()
    parser.add_argument('-f', dest='in_file_name', type=str, required=True)
    parser.add_argument('-y', dest='year_ranges', type=str, required=True)  # assume may overlap, (y1,y2 y3,y4 ...)
    parser.add_argument('-p', dest='position_ranges', type=str, required=True)  # assume to not overlap
    parser.add_argument('-o', dest='out_file_name', type=str, required=True)
    cmd_args = parser.parse_args()

    def parse_range(r):
        return [(int(i.split(',')[0]), int(i.split(',')[1])) for i in r.split(' ')]
    year_ranges = parse_range(cmd_args.year_ranges)
    pos_ranges = parse_range(cmd_args.position_ranges)

    with open(cmd_args.in_file_name, 'r') as f:
        r = reader(f)
        first_row = next(r)[POS_START_IND:]
        ind2pos = {}
        aa_counts = {}
        for i in range(0, len(first_row)):
            pos = int(first_row[i])
            if any(list(map(lambda y: y[0] <= pos <= y[1], pos_ranges))):
                ind2pos[i] = pos
                aa_counts[pos] = {aa: {r: 0 for r in year_ranges} for aa in AMINOACIDS}

        for row in r:


if __name__ == '__main__':
    main()


