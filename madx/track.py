import beamline

def read_madx_tracking(file):
    """Read a MAD-X Tracking onetable=true file to a dataframe."""
    column_names = ['ID','TURN','X','PX','Y','PY','T','PT','S','E']
    data = pd.read_csv(file, skiprows=54, delim_whitespace=True, names=column_names)
    return data.apply(pd.to_numeric, errors="ignore").dropna()


def read_ptc_tracking(file):
    """Read a PTC Tracking 'one' file to a dataframe."""
    column_names = ['ID', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E']
    data = pd.read_csv(file, skiprows=9, delim_whitespace=True,
                       names=column_names) \
              .apply(pd.to_numeric, errors="ignore").dropna()
    return data[data['TURN'] == 1]


def track(self, **kwargs):
    """Compute the distribution of the beam as it propagates through the beamline."""
    if kwargs.get('ptc', False):
        self.__flag_ptc = True
    m.beam()
    m.track(self.__beam.distribution, ptc=self.__flag_ptc)
    errors = m.run(self.context).fatals
    self.__madx_input = m.input
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
    if 'BEAM' in self.__beamline:
        self.__beamline.drop('BEAM', inplace=True, axis=1)
    self.__beamline = self.__beamline.merge(pd.DataFrame(tmp, columns=['BEAM']),
                                            left_on='AT_CENTER_TRUNCATED',
                                            right_index=True,
                                            how='left').sort_values(by='AT_CENTER')
    self.__beamline.drop('AT_CENTER_TRUNCATED', axis=1, inplace=True)
    self.__beamline.sort_values(by='S', inplace=True)
    return self.__beamline
