from argparse import ArgumentParser
from csv import reader, writer
from os.path import dirname, basename, join
from src.HIV.constants import AMINOACIDS


POS_START_IND = 4  # amino acid data starts at column 4
YEAR_IND = 1  # index of column holding year


def parse_cmdln_args():
    parser = ArgumentParser()
    parser.add_argument('-f', dest='in_file_name', type=str, required=True,
                        help='name of the input csv file. it should be a table with columns Country, Year, Patient, '
                             'Acession, position 1, position 2, ...')
    parser.add_argument('-y', dest='year_ranges', type=str, required=True,
                        help='year ranges to group results. it should be plain text with the format '
                             'range_1_start,range_1_end range_2_start,range_2_end ... end must be less than start, and '
                             'it is ok for the ranges to overlap')
    parser.add_argument('-p', dest='position_ranges', type=str, required=True,
                        help='positions to calculate the distribution for, in the same format as -y option, each '
                             'position will get a separate output file. These ranges should not overlap')
    parser.add_argument('-o', dest='out_file_name', type=str, required=True,
                        help='absolute path to output file. If "foo/bar/file.csv" is given, output file names will '
                             'look like foo/bar/position_1_file.csv, foo/bar/position_2_file.csv and so on')
    return parser.parse_args()


def parse_range(r):
    return [range(int(i.split(',')[0]), int(i.split(',')[1]) + 1) for i in r.split(' ')]


def setup(fn, year_ranges, pos_ranges):
    ind2pos = {}
    aa_counts = {}
    sums = {y: 0 for y in year_ranges}
    with open(fn, 'r') as f:
        r = reader(f)
        first_row = next(r)
        for i in range(POS_START_IND, len(first_row)):
            pos = int(first_row[i])
            if any(list(map(lambda x: pos in x, pos_ranges))):
                ind2pos[i] = pos
                aa_counts[pos] = {aa: {r: 0 for r in year_ranges} for aa in AMINOACIDS}

        for row in r:
            year = int(row[YEAR_IND])
            for year_range in year_ranges:
                if year in year_range:
                    for i in ind2pos:
                        if row[i] == '-':
                            continue
                        aa_counts[ind2pos[i]][row[i]][year_range] += 1
                    sums[year_range] += 1

    return aa_counts, sums


def write_results(ofn, aa_counts, year_ranges, sums):
    header = [''] + [f'{r.start}-{r.stop - 1}' for r in year_ranges]
    for pos in aa_counts:  # each position will get a separate file
        fn = join(dirname(ofn), str(pos) + '_' + basename(ofn))
        with open(fn, 'w') as f:
            w = writer(f, lineterminator='\n')
            w.writerow(header)
            for aa in aa_counts[pos]:
                row = [aa]
                for year_range in year_ranges:
                    row.append(str(aa_counts[pos][aa][year_range] / sums[year_range] * 100))
                w.writerow(row)


def main():
    cmd_args = parse_cmdln_args()
    year_ranges = parse_range(cmd_args.year_ranges)
    pos_ranges = parse_range(cmd_args.position_ranges)
    aa_counts, sums = setup(cmd_args.in_file_name, year_ranges, pos_ranges)
    write_results(cmd_args.out_file_name, aa_counts, year_ranges, sums)


if __name__ == '__main__':
    main()
