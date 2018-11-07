"""
Author: Rentian Dong

This script computes the distribution of amino acid profiles
"""


from argparse import ArgumentParser
from csv import reader, writer
from os.path import dirname, basename, join
from src.HIV.constants import AMINOACIDS


POS_START_IND = 1  # amino acid data starts at column 1
YEAR_IND = 0  # index of column holding year


def __parse_cmdln_args():
    """
    parse the cmdln arguments, use -h option for more help
    :return: an argparse.Namespace object with the parsed arguments
    """
    parser = ArgumentParser()
    parser.add_argument('-f', dest='in_file_name', type=str, required=True,
                        help='name of the input csv file. it should be a table with columns Year, position_1, '
                             'position _2, ...')
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


def __parse_range(ranges):
    """
    parse the cmdln input of -y option
    :param ranges: a string of the form 'y1,y2 y3,y4 ...'
    :return: [range(y1, y2), range(y3, y4) ...]
    """
    parsed_ranges = []
    for year_pair in ranges.split():
        [start, end] = year_pair.split(',')
        parsed_ranges.append(range(int(start), int(end) + 1))
    return parsed_ranges


def __map_ind2pos(file_name, pos_ranges):
    """
    return a mapping of positions (specified in cmdln args) to column indices
    :param file_name: name of csv file to create mapping from
    :param pos_ranges: a list of range objects, position contained in any of these range objects will be mapped
    :return: a dictionary of the form {index: position}
    """
    ind2pos = {}
    with open(file_name, 'r') as f:
        r = reader(f)
        first_row = next(r)
    for i in range(POS_START_IND, len(first_row)):
        pos = int(first_row[i])
        if any(map(lambda x: pos in x, pos_ranges)):  # include only positions supplied through cmdln args
            ind2pos[i] = pos
    return ind2pos


def __calc_distribution(file_name, year_ranges, ind2pos):
    """
    calculate the amino acid distribution within the given csv file
    :param file_name: name of the csv file, formatted as (Year, position_1, position _2, ...)
    :param year_ranges: a list of range objects, sequences whose sampling year is contained in any of these range
    objects will be included
    :param ind2pos: a dictionary of the form { column index: position }, obtained with map_ind2pos
    :return: a dictionary will percentage amino acid distributions. the returned dictionary has the format:
    { position: { amino acid: { year_range: percentage distribution } } }
    """

    # aa_counts: dictionary to keep track of frequency of each amino acid in each position and year range
    # aa_counts has the foramt { position: { amino acid: { year_range: frequency } } }
    aa_counts = {pos: {aa: {r: 0 for r in year_ranges} for aa in AMINOACIDS} for pos in ind2pos.values()}

    # sums: a dictionary to keep track of total
    sums = {year_range: 0 for year_range in year_ranges}

    with open(file_name) as file:
        r = reader(file)
        next(r)  # discard header row in file
        for row in r:
            year = int(row[YEAR_IND])
            for year_range in year_ranges:
                if year in year_range:
                    for i in ind2pos:
                        if row[i] == '-':
                            continue
                        aa_counts[ind2pos[i]][row[i]][year_range] += 1
                    sums[year_range] += 1

    # compute the amino acid distribution profile
    for pos in aa_counts:  # each position will get a separate file
        for aa in aa_counts[pos]:
            for year_range in year_ranges:
                aa_counts[pos][aa][year_range] = aa_counts[pos][aa][year_range] / sums[year_range] * 100

    return aa_counts


def __gen_out_file_name(pos, out_file_name):
    """
    construct an output file name
    :param pos: position associated with the distribution
    :param out_file_name: absolute path to output file name, used as a template.
    :return: the constructed file name
    """
    return join(dirname(out_file_name), str(pos) + '_' + basename(out_file_name))


def __build_header(aa_counts):
    """
    construct a header row to be written to output csv using the given distribution data
    :param aa_counts: distribution data returned from __calc_distribution
    :return: a list that looks like ['', 'year1-year2', 'year3-year4', ...]
    """
    year_ranges = next(iter(next(iter(aa_counts.values())).values())).keys()
    return [''] + [f'{r.start}-{r.stop - 1}' for r in year_ranges]


def __write_results(out_file_name, aa_counts):
    """
    write percentage distribution data to output file
    :param out_file_name: absolute path output file
    :param aa_counts: the dictionary containing percentage distribution data, obtained from calc_distribution
    :return: None, writes to the desginated file on disk
    """

    header = __build_header(aa_counts)
    for pos in aa_counts:  # each position will get a separate file
        with open(__gen_out_file_name(pos, out_file_name), 'w') as out_file:
            w = writer(out_file, lineterminator='\n')
            w.writerow(header)
            for aa in aa_counts[pos]:
                row = [aa]
                for year_range in aa_counts[pos][aa]:
                    row.append(str(aa_counts[pos][aa][year_range]))
                w.writerow(row)


def main():
    """
    the entry point of this script. the main method puts together all helper functions to calculate and then write
    amino acid distribution profile to an output file
    """
    cmd_args = __parse_cmdln_args()
    year_ranges = __parse_range(cmd_args.year_ranges)
    pos_ranges = __parse_range(cmd_args.position_ranges)
    ind2pos_mapping = __map_ind2pos(cmd_args.in_file_name, pos_ranges)
    distribution = __calc_distribution(cmd_args.in_file_name, year_ranges, ind2pos_mapping)
    __write_results(cmd_args.out_file_name, distribution)


if __name__ == '__main__':
    main()
