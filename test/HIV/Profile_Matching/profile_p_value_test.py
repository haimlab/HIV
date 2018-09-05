import unittest

import src.HIV.Profile_Matching.profile_p_value as p
from file_parse import get_all_static_profiles
from math import log10


class TestStringMethods(unittest.TestCase):

    def test_log_convert(self):
        self.assertEqual(p.logConvert(0), 0)
        self.assertEqual(p.logConvert(55), log10(55) + 2)


if __name__ == '__main__':
    unittest.main()
