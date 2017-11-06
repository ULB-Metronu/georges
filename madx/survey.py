import os
import pandas as pd
from .. import beamline
from .madx import Madx

MADX_SURVEY_HEADERS_SKIP_ROWS = 6
MADX_SURVEY_DATA_SKIP_ROWS = 8


class SurveyException(Exception):
    """Exception raised for errors in the Survey module."""

    def __init__(self, m):
        self.message = m


def read_survey(file):
    """Read a MAD-X survey file to a datraframe"""
    headers = pd.read_csv(file, skiprows=MADX_SURVEY_HEADERS_SKIP_ROWS, nrows=0, delim_whitespace=True)
    headers.drop(headers.columns[[0, 1]], inplace=True, axis=1)
    return pd.read_csv(file,
                       header=None,
                       names=headers.columns.values,
                       na_filter=False,
                       delim_whitespace=True,
                       skiprows=MADX_SURVEY_DATA_SKIP_ROWS
                       )


def survey(**kwargs):
    """Compute the survey of the beamline."""
    # Process arguments
    line = kwargs.get("line", None)
    if beamline is None:
        raise SurveyException("A beamline is expected.")
    m = Madx(beamlines=[line])
    m.beam(line.name)
    m.survey(start=kwargs.get('start', False))
    errors = m.run(**kwargs).fatals
    if kwargs.get("debug", False):
        print(m.input)
    if len(errors) > 0:
        print(errors)
        raise SurveyException("MAD-X ended with fatal error.")
    madx_survey = read_survey(os.path.join(".", 'survey.out'))
    line_with_survey = madx_survey.merge(line.line,
                                         left_index=True,
                                         right_index=True,
                                         how='outer',
                                         suffixes=('_SURVEY', '')
                                         ).sort_values(by='S')
    return beamline.Beamline(line_with_survey)
