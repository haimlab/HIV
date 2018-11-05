import unittest

import file_parse
from src.HIV import constants
from os.path import join
from math import log10
from copy import deepcopy

#TODO add a test for invalid filter property
import helpers


class TestFileParse(unittest.TestCase):

    def test_parse_file_name(self):

        # positive cases
        sol = 'C', 'EU', 1000
        self.assertEqual(file_parse.parse_file_name(join('data', 'static', 'C_EU_1000.csv')), sol)
        sol = 'AE', 'BRAZIL', 123
        self.assertEqual(file_parse.parse_file_name(join('data', 'static', 'AE_BRAZIL_123.csv')), sol)

        # negative cases
        self.assertRaises(ValueError, file_parse.parse_file_name, 'not-exist-dir')
        self.assertRaises(ValueError, file_parse.parse_file_name, 'not-exist-dir')

    def test_read_static(self):
        file_name = join('data', 'static', 'A1_ALL_295.csv')
        profile = file_parse.read_static(file_name)

        self.assertEqual(profile.clade, 'A1')
        self.assertEqual(profile.region, 'ALL')
        self.assertEqual(profile.position, 295)

    def test_get_all_static_profiles(self):
        all_p = file_parse.get_all_static_profiles()
        self.assertEqual(len(all_p.profiles), 125)

    def test_calc_year(self):
        self.assertEqual(helpers.calc_year('[2010, 2015]'), (2010 + 2015) / 2)
        self.assertEqual(helpers.calc_year('[2000, 2007]'), (2000 + 2007) / 2)

    def test_filter(self):
        all_p = file_parse.get_all_static_profiles()
        pos_295 = all_p.filter(295)
        pos_332 = all_p.filter(332)
        pos_339 = all_p.filter(339)
        pos_392 = all_p.filter(392)
        pos_448 = all_p.filter(448)
        for p in pos_332.profiles:
            self.assertTrue(p.position == 332)
        l_295 = len(pos_295.profiles)
        l_332 = len(pos_332.profiles)
        l_339 = len(pos_339.profiles)
        l_392 = len(pos_392.profiles)
        l_448 = len(pos_448.profiles)
        my_sum = l_295 + l_332 + l_339 + l_392 + l_448
        pos_2F5 = [662, 663, 664, 665, 667]
        for p in pos_2F5:
            my_sum += len(all_p.filter(p).profiles)
        l_tot = len(all_p.profiles)
        self.assertEqual(my_sum, l_tot)

    def test_filter_with_named_args(self):
        all_static_profs = file_parse.get_all_static_profiles()
        sub = all_static_profs.filter('B', 295)
        for s in sub.profiles:
            self.assertEqual(s.position, 295)
            self.assertEqual(s.clade, 'B')

    def test_filter_with_vargs(self):
        all_static_profs = file_parse.get_all_static_profiles()
        sub = all_static_profs.filter('B', 295)
        for s in sub.profiles:
            self.assertEqual(s.position, 295)
            self.assertEqual(s.clade, 'B')

    def test_filter_with_mixed_args(self):
        all_static_profs = file_parse.get_all_static_profiles()
        sub = all_static_profs.filter(295, 'B')
        for s in sub.profiles:
            self.assertEqual(s.position, 295)
            self.assertEqual(s.clade, 'B')

    def test_shuffle(self):
        all_p = file_parse.get_all_static_profiles()
        shuffled_p = deepcopy(all_p)
        shuffled_p.shuffle('clade')
        all_dict = {c: 0 for c in constants.CLADES}
        shuffled_dict = {c: 0 for c in constants.CLADES}
        for a, b in zip(all_p.profiles, shuffled_p.profiles):
            all_dict[a.clade] += 1
            shuffled_dict[b.clade] += 1
        self.assertEqual(all_dict, shuffled_dict)

    def test_log_convert(self):
        p = file_parse.StaticProfile(
            'A1',
            'ALL',
            100
        )
        p.distr['A'] =  20
        p.distr['C'] = 0
        p.distr['D'] = 100
        p = p.log_convert()
        self.assertAlmostEqual(2.301, p.distr['A'], delta=0.01)
        self.assertEqual(0, p.distr['C'])
        self.assertAlmostEqual(3, p.distr['D'], delta=0.01)

    def test_log_convert_helper(self):
        for i in range(11, 30):
            self.assertEqual(0, helpers.log_convert(1 / i))
        for i in range(10, 100):
            self.assertAlmostEqual(log10(i) + 1, helpers.log_convert(i), delta=0.01)

    def test_attr_list(self):

        # dynamic profile lists
        p = file_parse.get_all_dynamic_profiles()
        attrs = set(p.attr_list('clade'))
        self.assertEqual(set(constants.CLADES), attrs)

        # static profile lists
        p = file_parse.get_all_static_profiles()
        p = p.filter('EU', 'B')
        self.assertEqual(['EU'], p.attr_list('region'))
        self.assertCountEqual(['B'], p.attr_list('clade'))
