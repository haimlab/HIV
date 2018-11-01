from argparse import ArgumentParser

AMINOACIDS = ('Z', 'N', 'T', 'S', 'D', 'E',
              'K', 'R', 'H', 'Y', 'Q', 'I',
              'L', 'V', 'A', 'C', 'F', 'G',
              'M', 'P', 'W')

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
    positions = parse_range(cmd_args.position_ranges)


if __name__ == '__main__':
    main()

    # aa_counts = {p: {r: {aa: 0 for aa in AMINOACIDS} for r in year_ranges} for p in positions}
    #
    # with open(ifn, 'r') as f:
    #     fr = next(f)
    #     positions = filter(lambda x: any(list(map(lambda y: y[0] <= int(x) <= y[1], year_ranges))), fr)
    #     print(positions)
