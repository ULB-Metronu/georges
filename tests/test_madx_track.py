import unittest

from georges import madx


class TestMadxTrack(unittest.TestCase):

    def test_invalid_instanciation(self):
        with self.assertRaisesRegex(georges.madx.tracking.TrackException,
                                    "Beamline, Beam and MAD-X objects need to be defined."):
            madx.track()
