import unittest

import file_parse
import constants

class TesTFileParse(unittest.TestCase):

    def test_parse_file_name(self):

        # positive cases
        sol = constants.Clade.C, constants.Region.EU, 1000
        self.assertEqual(file_parse.parse_file_name('data\\static\\C_EU_1000.csv'), sol)
        sol = constants.Clade.AE, constants.Region.BRAZIL, 123
        self.assertEqual(file_parse.parse_file_name('data\\static\\AE_BRAZIL_123.csv'), sol)

        # negative cases
        self.assertRaises(ValueError, file_parse.parse_file_name, 'balh\\JACK_BRAZIL_234.csv')
        self.assertRaises(ValueError, file_parse.parse_file_name, 'blah\\AE_THISISSPARTA!_123.csv')

    def test_read_static(self):
        file_name = 'data\\static\\A1_ALL_295.csv'
        profile = file_parse.read_static(file_name)

        self.assertEqual(profile.clade(), constants.Clade.A1)
        self.assertEqual(profile.region(), constants.Region.ALL)
        self.assertEqual(profile.position(), 295)

    def test_get_all_static_profiles(self):
        all_p = file_parse.get_all_static_profiles()
        self.assertEqual(len(all_p.get_all_profiles()), 75)

    def test_calc_year(self):
        self.assertEqual(file_parse.calcYear('[2010, 2015]'), (2010 + 2015) / 2)
        self.assertEqual(file_parse.calcYear('[2000, 2007]'), (2000 + 2007) / 2)

    def test_filter(self):
        all_p = file_parse.get_all_static_profiles()
        pos_332 = all_p.filter(position=332)
        for p in pos_332.get_all_profiles():
            self.assertTrue(p.position() == 332)
        for p in all_p.get_all_profiles():
            if p not in pos_332.get_all_profiles():
                self.assertTrue(p.position() != 332)
