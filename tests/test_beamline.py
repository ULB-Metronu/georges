import unittest

import pandas as pd

from georges import beamline

PREFIX = 'tests/test_data'


class TestBeamline(unittest.TestCase):

    def test_invalid_instanciation_empty(self):
        with self.assertRaisesRegex(beamline.BeamlineException, 'Single argument expected.'):
            beamline.Beamline()

    def test_invalid_instanciation_no_beamline(self):
        with self.assertRaisesRegex(OSError, 'File .* does not exist'):
            beamline.Beamline('SOME_BEAMLINE_NOT_IN_PATH')

    def test_invalid_instanciation_multiple_file_arguments(self):
        with self.assertRaisesRegex(beamline.BeamlineException, "Single argument expected."):
            beamline.Beamline("FILE1", "FILE2")

    def test_invalid_instanciation_multiple_dataframes_arguments(self):
        with self.assertRaises(beamline.BeamlineException):
            beamline.Beamline(pd.DataFrame(), pd.DataFrame())

    def test_invalid_instanciation_empty_file(self):
        with self.assertRaisesRegex(pd.io.common.EmptyDataError, 'No columns to parse from file'):
            beamline.Beamline('empty', prefix=PREFIX)

    def test_invalid_instanciation_empty_dataframe(self):
        with self.assertRaisesRegex(beamline.BeamlineException, "Empty dataframe."):
            beamline.Beamline(pd.DataFrame())

    def test_valid_instanciation_file(self):
        b = beamline.Beamline('test_beamline', prefix=PREFIX)
        self.assertEqual(b.length, 0.215)

    def test_valid_survey_conversion_trigger(self):
        self.assertIs(beamline.Beamline('test_beamline', prefix=PREFIX).converted_from_survey, False)
        self.assertIs(beamline.Beamline('test_beamline_survey', prefix=PREFIX).converted_from_survey, True)

    def test_invalid_survey_data(self):
        with self.assertRaisesRegex(beamline.BeamlineException, "Trying to infer sequence from survey data: X and Y must be provided."):
            beamline.Beamline('test_beamline_invalid_survey', prefix=PREFIX)

    def test_invalid_elements_data_empty_list(self):
        with self.assertRaisesRegex(beamline.BeamlineException, "Invalid data type for 'elements'"):
            beamline.Beamline('test_beamline', elements=[], prefix=PREFIX)

    def test_invalid_elements_data_file(self):
        with self.assertRaisesRegex(OSError, "File .* does not exist"):
            beamline.Beamline('test_beamline', elements="nonexisting_file")

    def test_name(self):
        self.assertEqual(beamline.Beamline('test_beamline', prefix=PREFIX).name, 'TEST_BEAMLINE')

    def test_length(self):
        self.assertEqual(beamline.Beamline('test_beamline', prefix=PREFIX).length, 0.215)

if __name__ == '__main__':
    unittest.main()
