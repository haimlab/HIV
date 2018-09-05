import unittest

import src.HIV.Profile_Matching.profile_p_value as p
from file_parse import get_all_static_profiles


class TesTFileParse(unittest.TestCase):

    def test_read_static(self):
        file_name =