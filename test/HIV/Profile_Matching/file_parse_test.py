import unittest

import file_parse
import constants


class TestFileParse(unittest.TestCase):

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
        self.assertEqual(len(all_p.get_all_profiles()), 125)

    def test_calc_year(self):
        self.assertEqual(file_parse.calcYear('[2010, 2015]'), (2010 + 2015) / 2)
        self.assertEqual(file_parse.calcYear('[2000, 2007]'), (2000 + 2007) / 2)

    def test_filter(self):
        all_p = file_parse.get_all_static_profiles()
        pos_295 = all_p.filter(295)
        pos_332 = all_p.filter(332)
        pos_339 = all_p.filter(339)
        pos_392 = all_p.filter(392)
        pos_448 = all_p.filter(448)
        for p in pos_332.get_all_profiles():
            self.assertTrue(p.position() == 332)
        l_295 = len(pos_295.get_all_profiles())
        l_332 = len(pos_332.get_all_profiles())
        l_339 = len(pos_339.get_all_profiles())
        l_392 = len(pos_392.get_all_profiles())
        l_448 = len(pos_448.get_all_profiles())
        my_sum = l_295 + l_332 + l_339 + l_392 + l_448
        pos_2F5 = [662, 663, 664, 665, 667]
        for p in pos_2F5:
            my_sum += len(all_p.filter(p).get_all_profiles())
        l_tot = len(all_p.get_all_profiles())
        self.assertEqual(my_sum, l_tot)

    def test_filter_with_named_args(self):
        all_static_profs = file_parse.get_all_static_profiles()
        sub = all_static_profs.filter(constants.Clade.B, 295)
        for s in sub.get_all_profiles():
            self.assertEqual(s.position(), 295)
            self.assertEqual(s.clade(), constants.Clade.B)

    def test_filter_with_vargs(self):
        all_static_profs = file_parse.get_all_static_profiles()
        sub = all_static_profs.filter(constants.Clade.B, 295)
        for s in sub.get_all_profiles():
            self.assertEqual(s.position(), 295)
            self.assertEqual(s.clade(), constants.Clade.B)

    def test_filter_with_mixed_args(self):
        all_static_profs = file_parse.get_all_static_profiles()
        sub = all_static_profs.filter(295, constants.Clade.B)
        for s in sub.get_all_profiles():
            self.assertEqual(s.position(), 295)
            self.assertEqual(s.clade(), constants.Clade.B)

    def test_shuffle(self):
        all_p = file_parse.get_all_static_profiles()
        shuffled_p = all_p.shuffle(constants.FilterProperties.CLADE)
        all_dict = {c: 0 for c in constants.Clade}
        shuffled_dict = {c: 0 for c in constants.Clade}
        for a, b in zip(all_p.get_all_profiles(), shuffled_p.get_all_profiles()):
            all_dict[a.clade()] += 1
            shuffled_dict[b.clade()] += 1
        self.assertEqual(all_dict, shuffled_dict)



