import unittest
import pandas as pd
from context import beamline


class TestBeamline(unittest.TestCase):

    def test_invalid_instanciation_empty(self):
        with self.assertRaisesRegex(beamline.BeamlineException, 'No beamline defined.'):
            beamline.Beamline()

    def test_invalid_instanciation_no_beamline(self):
        with self.assertRaisesRegex(OSError, 'File .* does not exist'):
            beamline.Beamline('SOME_BEAMLINE_NOT_IN_PATH')

    def test_invalid_instanciation_multiple_dataframes(self):
        with self.assertRaises(beamline.BeamlineException):
            beamline.Beamline(pd.DataFrame(), pd.DataFrame())

    def test_invalid_instanciation_empty_file(self):
        with self.assertRaisesRegex(pd.io.common.EmptyDataError, 'No columns to parse from file'):
            beamline.Beamline('empty')

    def test_invalid_instanciation_empty_dataframe(self):
        with self.assertRaisesRegex(beamline.BeamlineException, "Empty dataframe."):
            beamline.Beamline(pd.DataFrame())

    def test_valid_instanciation_file(self):
        b = beamline.Beamline('test_beamline')
        self.assertEqual(b.length, 0.215)

    def test_valid_survey_conversion_trigger(self):
        self.assertIs(beamline.Beamline('test_beamline').converted_from_survey, False)
        self.assertIs(beamline.Beamline('test_beamline_survey').converted_from_survey, True)

    def test_invalid_survey_data(self):
        with self.assertRaisesRegex(beamline.BeamlineException, "Trying to infer sequence from survey data: X and Y must be provided."):
                beamline.Beamline('test_beamline_invalid_survey')

    def test_name(self):
        self.assertEqual(beamline.Beamline('test_beamline').name, 'TEST_BEAMLINE')

    def test_length(self):
        self.assertEqual(beamline.Beamline('test_beamline').length, 0.215)

if __name__ == '__main__':
    unittest.main()