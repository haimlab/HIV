import unittest
import file_parse
import constants
import profile_p_value


class TestStringMethods(unittest.TestCase):

    def test_parse_file_name(self):

        all_prof = file_parse.get_all_static_profiles()
        sub_prof = profile_p_value.select_sub_group(all_prof, ['B,NA', 'C,EU', 'AE,TH'], [295, 667])
        for p in sub_prof.get_all_profiles():
            self.assertTrue(p.position() == 295 or p.position() == 667)
        for p in sub_prof.get_all_profiles():
            self.assertTrue(
                p.clade() == constants.Clade.B and p.region() == constants.Region.NA or
                p.clade() == constants.Clade.C and p.region() == constants.Region.EU or
                p.clade() == constants.Clade.AE and p.region() == constants.Region.TH
            )