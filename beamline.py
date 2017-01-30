import sys
import os.path
import pandas as pd
import numpy as np
import numpy.linalg as npl
import georges.beam as beam
import georges.madx.madx as madx

DEFAULT_EXT = 'csv'


def compute_derived_data(row):
    """Compute derived data for a beamline element."""
    # Add a bunch of mandatory columns
    for c in ['E1', 'E2', 'FINT', 'ANGLE', 'HGAP', 'THICK', 'TILT']:
        if pd.isnull(row.get(c)):
            row[c] = np.nan

    # Corner case
    if pd.isnull(row.get('CLASS')):
        row['CLASS'] = row.get('TYPE', 'MARKER')

    # Should be changed to infer the list of the module functions
    for d in row.index | ['AT_ENTRY', 'AT_EXIT', 'AT_CENTER', 'ORBIT_LENGTH', 'ANGLE', 'E1', 'E2']:
        if pd.isnull(row.get(d)) and not pd.isnull(row.get(d + '_ELEMENT')):
            row[d] = row[d + '_ELEMENT']
        if pd.isnull(row.get(d)):
            row[d] = getattr(sys.modules[__name__], d.lower(), lambda r: np.nan)(row)

    return row


def at_entry(r):
    """Try to compute the element's entry 's' position from other data."""
    if not pd.isnull(r.get('AT_CENTER')):
        return r['AT_CENTER'] - r['ORBIT_LENGTH']/2.0
    elif not pd.isnull(r.get('AT_EXIT')):
        return r['AT_EXIT'] - r['ORBIT_LENGTH']
    else:
        return np.nan


def at_center(r):
    """Try to compute the element's center 's' position from other data."""
    if pd.isnull(r.get('ORBIT_LENGTH')):
        return np.nan
    if not pd.isnull(r.get('AT_ENTRY')):
        return r['AT_ENTRY'] + r['ORBIT_LENGTH']/2.0
    elif not pd.isnull(r.get('AT_EXIT')):
        return r['AT_EXIT'] - r['ORBIT_LENGTH']/2.0
    else:
        return np.nan


def at_exit(r):
    """Try to compute the element's entry 's' position from other data."""
    if pd.isnull(r.get('ORBIT_LENGTH')):
        return np.nan
    if not pd.isnull(r.get('AT_ENTRY')):
        return r['AT_ENTRY'] + r['ORBIT_LENGTH']
    elif not pd.isnull(r.get('AT_CENTER')):
        return r['AT_CENTER'] + r['ORBIT_LENGTH']/2.0
    else:
        return np.nan


def orbit_length(r):
    """Try to compute the element's orbit length from other data."""
    if pd.isnull(r.get('LENGTH')):
        return 0.0
    if pd.isnull(r.get('ANGLE')):
        return r['LENGTH']
    return r['ANGLE']*np.pi/180.0 * r['LENGTH'] / (2.0 * np.sin((r['ANGLE']*np.pi/180.0) / 2.0))


class BeamlineException(Exception):
    """Exception raised for errors in the Beamline module."""

    def __init__(self, m):
        self.message = m


class Beamline:
    """A beamline or accelerator model.

    The internal representation is essentially a set of pandas DataFrame.
    """

    def __init__(self, *args, **kwargs):
        # We need the kwargs first to parse the args correctly
        self.__path = kwargs.get('path', '.')
        self.__prefix = kwargs.get('prefix', '')
        self.__elements = kwargs.get('elements', None)
        self.__length = kwargs.get('length', 0)
        self.__strengths = None
        self.__beam = kwargs.get('beam', None)
        self.__flag_ptc = kwargs.get('ptc', False)

        # Some type inference to get the elements right
        # Elements as a file name
        if self.__elements and isinstance(self.__elements, str):
            self.__build_elements_from_file(self.__elements)
        # Elements as a list to be converted onto a DataFrame
        elif self.__elements and isinstance(self.__elements, list):
            self.__elements = pd.DataFrame(self.__elements)

        # Now let's process the args
        for i, arg in enumerate(args):
            # Some type inference to get the sequence right
            # Sequence from file name
            if isinstance(arg, str) and i == 0:
                self.__name = arg.upper()
                self.__build_from_files([arg])
            # Sequence from a list of files
            if isinstance(arg, list) and len(arg) > 0 and i == 0:
                self.__name = '_'.join(arg).upper()
                self.__build_from_files(arg)
            # Sequence from a pandas.DataFrame
            if isinstance(arg, pd.DataFrame) and i == 0:
                self.__name = getattr(arg, 'name', 'BEAMLINE')
                self.__beamline = arg

        # Check before hand if the survey will need to be converted
        survey = False
        if self.__beamline.get(['AT_ENTRY', 'AT_CENTER', 'AT_EXIT']) is None:
            if self.__beamline.get(['X', 'Y']) is not None:
                survey = True

        # Expand elements onto MAD-X native elements
        if self.__elements is not None and self.__beamline is not None:
            self.__expand_elements_data()

        # Compute derived data until a fixed point sequence is reached
        self.__expand_sequence_data()

        # Infer if the sequence is given as a survey and convert to positions
        if survey:
            self.__convert_survey_to_sequence()
            # Re-expand
            self.__expand_sequence_data()

        # Compute the sequence length if needed
        if self.__length == 0:
            self.__length = self.__beamline['AT_EXIT'].values[-1]

        # Flag to distinguish MAD-X generated elements from beamline elements
        self.__beamline['PHYSICAL'] = True


    @property
    def name(self):
        """The sequence name."""
        return self.__name

    @property
    def length(self):
        """The sequence length."""
        return self.__length

    @property
    def strengths(self):
        """Strengths of the various elements of the beamline."""
        return self.__strengths

    @strengths.setter
    def strengths(self, strengths):
        if not strengths.index == self.__strengths:
            raise BeamlineException("Trying to set the strengths for an invalid elements list.")
        self.__strengths = strengths

    @property
    def elements(self):
        """Elements composing the beamline."""
        return self.__elements

    @property
    def line(self):
        """The beamline representation."""
        self.__beamline.name = self.name
        self.__beamline.length = self.length
        return self.__beamline

    @property
    def twiss(self, **kwargs):
        """Compute the Twiss parameters of the beamline."""
        if kwargs.get('ptc', False):
            self.__flag_ptc = True
        m = madx.Madx(beamline=self, path=self.__path, madx='/usr/local/bin/madx-dev')
        m.beam()
        m.twiss(ptc=self.__flag_ptc)
        errors = m.run(self.__get_context()).fatals
        if len(errors) > 0:
            print(m.input)
            print(errors)
            raise BeamlineException("MAD-X ended with fatal error.")
        madx_twiss = madx.read_madx_twiss(os.path.join(self.__path, 'twiss.outx'))
        self.__beamline = madx_twiss.merge(self.__beamline,
                                           left_index=True,
                                           right_index=True,
                                           how='outer',
                                           suffixes=('_TWISS', '')
                                           ).sort_values(by='S')
        return self.__beamline

    @property
    def track(self, **kwargs):
        """Compute the distribution of the beam as it propagates through the beamline."""
        if kwargs.get('ptc', False):
            self.__flag_ptc = True
        m = madx.Madx(beamline=self, path=self.__path, madx='/usr/local/bin/madx-dev')
        m.beam()
        m.track(self.__beam.distribution, ptc=self.__flag_ptc)
        errors = m.run(self.__get_context()).fatals
        if len(errors) > 0:
            print(m.input)
            print(errors)
            raise BeamlineException("MAD-X ended with fatal error.")
        if self.__flag_ptc:
            madx_track = madx.read_ptc_tracking(os.path.join(self.__path, ''))
        else:
            madx_track = madx.read_madx_tracking(os.path.join(self.__path, 'tracking.outxone')).dropna()
            madx_track['PY'] = pd.to_numeric(madx_track['PY'])
        madx_track['S'] = round(madx_track['S'], 8)
        tmp = madx_track.query('TURN == 1').groupby('S').apply(lambda g: beam.Beam(g[['X', 'PX', 'Y', 'PY', 'PT']]))
        self.__beamline['AT_CENTER_TRUNCATED'] = round(self.__beamline['AT_CENTER'], 8)
        self.__beamline = self.__beamline.merge(pd.DataFrame(tmp, columns=['BEAM']),
                                                left_on='AT_CENTER_TRUNCATED',
                                                right_index=True,
                                                how='left').sort_values(by='S')
        self.__beamline.drop('AT_CENTER_TRUNCATED', axis=1, inplace=True)
        return self.__beamline

    def __get_context(self):
        return {
            'PARTICLE': self.__beam.particle,
            'PC': self.__beam.pc,
            'BETAREL': self.__beam.betarel,
            'BETAX': 0.0846155,
            'BETAY': 0.0846155,
            'ALPHAX': 0.0,
            'ALPHAY': 0.0,
            'DELTAP': 0.0,
            'DPP': 0.5e-2,
            'N_TRACKING': 5000,
            'EMITX': 14.3e-6,
            'EMITY': 14.3e-6,
        }

    def __build_from_files(self, names):
        files = [os.path.splitext(n)[0] + '.' + (os.path.splitext(n)[1] or DEFAULT_EXT) for n in names]
        sequences = [
            pd.read_csv(os.path.join(self.__path, self.__prefix, f), index_col = 'NAME') for f in files
            ]
        self.__beamline = pd.concat(sequences)

    def __expand_sequence_data(self):
        """Apply sequence transformation until a fixed point is reached."""
        tmp = self.__beamline.apply(compute_derived_data, axis=1)
        tmp2 = tmp
        while True:
            tmp, tmp2 = tmp2, tmp.apply(compute_derived_data, axis=1)
            if tmp.equals(tmp2):
                break
        self.__beamline = tmp2

    def __build_elements_from_file(self, file):
        file = os.path.splitext(file)[0] + '.' + (os.path.splitext(file)[1] or DEFAULT_EXT)
        self.__elements = pd.read_csv(os.path.join(self.__path, file), index_col='NAME')

    def __expand_elements_data(self):
        if self.__elements is None:
            return
        self.__beamline = self.__beamline.merge(self.__elements,
                                                left_on='TYPE',
                                                right_index=True,
                                                how='left',
                                                suffixes=('', '_ELEMENT')
                                                )
        # Angle conversion
        self.__beamline['ANGLE'] = self.__beamline['ANGLE'] / 180.0 * np.pi

    def __convert_survey_to_sequence(self):
        s = self.__beamline
        offset = s['ORBIT_LENGTH'][0] / 2.0
        if pd.isnull(offset):
            offset = 0
        self.__beamline['AT_CENTER'] = pd.DataFrame(npl.norm([
            s['X'].diff().fillna(0.0),
            s['Y'].diff().fillna(0.0)
        ], axis=0) - (
            s['LENGTH'].fillna(0.0) / 2.0 - s['ORBIT_LENGTH'].fillna(0.0) / 2.0
        ) + (
            s['LENGTH'].shift(1).fillna(0.0) / 2.0 - s['ORBIT_LENGTH'].shift(1).fillna(0.0) / 2.0
        )).cumsum() / 1000.0 + offset
