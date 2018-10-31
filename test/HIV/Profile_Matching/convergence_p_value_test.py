import unittest
from constants import AMINOACIDS
import helpers


class TestProfilePValue(unittest.TestCase):

    def test_envlopes_to_profile(self):
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

    def test_envelopes_to_profile(self):
        envelopes = []