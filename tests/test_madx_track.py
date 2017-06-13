import unittest
import madx.tracking


class TestMadxTrack(unittest.TestCase):

    def test_invalid_instanciation(self):
        with self.assertRaisesRegex(madx.tracking.TrackException,
                                    "Beamline, Beam and MAD-X objects need to be defined."):
            madx.track()
