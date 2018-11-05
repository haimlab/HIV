import unittest
import file_parse
import profile_p_value


class TestProfilePValue(unittest.TestCase):

    def test_select_sub_group(self):

        all_prof = file_parse.get_all_static_profiles()
        sub_prof = profile_p_value.select_sub_group(all_prof, ['B,NA', 'C,EU', 'AE,TH'], [295, 667])
        for p in sub_prof.profiles:
            self.assertTrue(p.position == 295 or p.position == 667)
        for p in sub_prof.profiles:
            self.assertTrue(
                p.clade == 'B' and p.region == 'NA' or
                p.clade == 'C' and p.region == 'EU' or
                p.clade == 'AE' and p.region == 'TH'
            )

    def test_select_sub_group_exclusions(self):

        all_prof = file_parse.get_all_static_profiles()
        sub_prof = profile_p_value.select_sub_group(all_prof, ['C,EU', 'AE,TH'], [295, 332, 392])
        num_c = 0
        num_ae = 0
        for p in sub_prof.profiles:
            if p.clade == 'C':
                self.assertTrue(p.position != 295)
                num_c += 1
            if p.clade == 'AE':
                self.assertTrue(p.position != 332)
                num_ae += 1
        self.assertEqual(2, num_ae)
        self.assertEqual(2, num_c)
