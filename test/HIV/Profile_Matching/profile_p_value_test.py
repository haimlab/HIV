import unittest
import file_parse
import constants


class TestStringMethods(unittest.TestCase):

    def test_filter_with_named_args(self):
        all_static_profs = file_parse.get_all_static_profiles()
        sub = all_static_profs.filter(clade=constants.Clade.B, position=295)
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
        sub = all_static_profs.filter(295, clade=constants.Clade.B)
        for s in sub.get_all_profiles():
            self.assertEqual(s.position(), 295)
            self.assertEqual(s.clade(), constants.Clade.B)
