import os, re, io
import pandas as pd
import numpy as np
import ROOT
import root_numpy
from .. import beamline
from .. import beam
from .bdsim import BDSim
from .. import physics


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def convert_array(data_frame, element):  # A mettre sous une meilleure forme

    new_df = pd.DataFrame()
    new_idx = 0
    for index, line in data_frame.iterrows():

        if np.size(line[element + ".x"]) > 0:
            new_df.at[new_idx, "X"] = line[element + ".x"][0]
            new_df.at[new_idx, "Y"] = line[element + ".y"][0]
            new_df.at[new_idx, "PX"] = line[element + ".xp"][0]
            new_df.at[new_idx, "PY"] = line[element + ".yp"][0]
            new_df.at[new_idx, "ENERGY"] = 1000 * (line[element + ".energy"][0] - 0.938)
            new_df.at[new_idx, "PARENTID"] = line[element + ".parentID"][0]
            new_idx += 1

    new_df['P'] = physics.energy_to_momentum(new_df['ENERGY'])
    return new_df


def read_tracking(element, evttree):
    """Read a BDSIM Tracking 'one' file to a dataframe."""

    if element['TYPE'] not in ['MARKER', 'SOLIDS']:
        data_element = root_numpy.tree2array(evttree, branches=[element.name + ".x",
                                                                element.name + ".y",
                                                                element.name + ".xp",
                                                                element.name + ".yp",
                                                                element.name + ".energy",
                                                                element.name + ".parentID"])
        df = pd.DataFrame(data=data_element)
        corrected_df = convert_array(df, element.name)
        tmp = beam.Beam(corrected_df[['X', 'PX', 'Y', 'PY', 'P']])
        return tmp


def track(**kwargs):
    """Compute the distribution of the beam as it propagates through the beamline..
    :param kwargs: parameters are:
        - line: the beamline on which twiss will be run
        - context: the associated context on which MAD-X is run
    """
    # Process arguments
    line = kwargs.get('line', None)
    b = kwargs.get('beam', None)
    context = kwargs.get('context', None)
    options = kwargs.get('options', None)

    if line is None or b is None or context is None:
        raise TrackException("Beamline, Beam, Context and BDsim objects need to be defined.")

    if options is None:
        print("Warning : no options is provided")

    bd = BDSim(beamlines=[line], **kwargs)

    # Write the input file for bdsim
    bd_beam = b.distribution.copy()
    p0 = physics.energy_to_momentum(b.energy)
    bd.track(bd_beam, p0)

    # Add options for Bdsim
    bd.set_options(options)

    # Run bdsim
    errors = bd.run(**kwargs).fatals

    # Create a new beamline to include the results
    l = line.line.copy()

    # Open the ROOT file and get the events
    f = ROOT.TFile('output.root')
    evttree = f.Get("Event")
    #
    # # Add columns which contains datas
    l['BEAM'] = l.apply(lambda g: read_tracking(g, evttree), axis=1)
    #
    return beamline.Beamline(l)
