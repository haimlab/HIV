import unittest
from src.HIV.constants import AMINOACIDS
import helpers


class TestProfilePValue(unittest.TestCase):

    def test_envlopes_to_profile_log(self):

        sequences = \
            ['Z'] * 10 + \
            ['N'] * 50 + \
            ['T'] * 1
        prof = helpers.envelopes_to_profile(sequences, 0)
        for aa in AMINOACIDS:
            if aa == 'Z':
                self.assertAlmostEqual(0.21467016498, prof[aa])
            elif aa == 'N':
                self.assertAlmostEqual(0.91364016932, prof[aa])
            else:
                self.assertEqual(0, prof[aa])

    def test_envelopes_to_profile_non_log(self):
        sequences = \
            ['Z'] * 10 + \
            ['N'] * 50 + \
            ['T'] * 1
        prof = helpers.envelopes_to_profile(sequences, 0, log=False)
        for aa in AMINOACIDS:
            if aa == 'Z':
                self.assertAlmostEqual(0.163934426229508196721311475409836065573770491803278688524, prof[aa])
            elif aa == 'N':
                self.assertAlmostEqual(0.819672131147540983606557377049180327868852459016393442622, prof[aa])
            elif aa == 'T':
                self.assertAlmostEqual(0.016393442622950819672131147540983606557377049180327868852, prof[aa])
            else:
                self.assertEqual(0, prof[aa])