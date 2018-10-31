import unittest
import file_parse
import constants
import profile_p_value


class TestProfilePValue(unittest.TestCase):

    def test_select_sub_group(self):

        all_prof = file_parse.get_all_static_profiles()
        sub_prof = profile_p_value.select_sub_group(all_prof, ['B,NA', 'C,EU', 'AE,TH'], [295, 667])
        for p in sub_prof.get_all_profiles():
            self.assertTrue(p.position == 295 or p.position == 667)
        for p in sub_prof.get_all_profiles():
            self.assertTrue(
                p.clade == constants.Clade.B and p.region == constants.Region.NA or
                p.clade == constants.Clade.C and p.region == constants.Region.EU or
                p.clade == constants.Clade.AE and p.region == constants.Region.TH
            )

    def test_select_sub_group_exclusions(self):

        all_prof = file_parse.get_all_static_profiles()
        sub_prof = profile_p_value.select_sub_group(all_prof, ['C,EU', 'AE,TH'], [295, 332, 392])
        num_c = 0
        num_ae = 0
        for p in sub_prof.get_all_profiles():
            if p.clade == constants.Clade.C:
                self.assertTrue(p.position != 295)
                num_c += 1
            if p.clade == constants.Clade.AE:
                self.assertTrue(p.position != 332)
                num_ae += 1
        self.assertEqual(num_ae, 2)
        self.assertEqual(num_c, 2)
