from argparse import ArgumentParser


# take a 1D list of years and put them into groups of two
def format_period_ranges(years):
    if len(years) % 2 != 0:
        raise Exception('cannot group odd length lists')
    return [(years[i], years[i + 1]) for i in range(len(years) - 1)]


