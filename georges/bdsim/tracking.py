import os, re, io
import pandas as pd
import numpy as np
import ROOT
import root_numpy
from .. import beamline
from .bdsim import BDSim
from .. import physics
from . import beam_bdsim


class TrackException(Exception):
    """Exception raised for errors in the Track module."""

    def __init__(self, m):
        self.message = m


def get_array_value(x):

    # TODO : optimize using the function from records of dataframe.

    if len(x) > 0:
        return x[0]
    else:
        return None


def read_tracking(element, evttree, **kwargs):
    """Read a BDSIM Tracking 'one' file to a dataframe."""

    if element['TYPE'] not in ['MARKER', 'SOLIDS', 'DEGRADER']:
        data_element = root_numpy.tree2array(evttree, branches=[element.name + ".x",
                                                                element.name + ".y",
                                                                element.name + ".xp",
                                                                element.name + ".yp",
                                                                element.name + ".energy",
                                                                element.name + ".parentID",
                                                                element.name + ".partID",
                                                                element.name + ".weight"])
        df = pd.DataFrame(data=data_element)
        df = df.applymap(lambda x: get_array_value(x))

        # TODO Improve the selection for the transport change for s position but wait for Laurie
        if kwargs.get("only_transport", False):
            exit = element['AT_EXIT']
            df['part_number'] = list(range(0, len(df)))
            prim = root_numpy.tree2array(evttree, branches=["PrimaryFirstHit.Z"])
            df1 = pd.DataFrame(data=prim)
            df1 = df1.applymap(lambda x: get_array_value(x))
            df1["part_number"] = list(range(0, len(df1)))
            df1["pos_z"] = df1["PrimaryFirstHit.Z"]
            # Only get particles which suffer an interaction before the element
            df1.query("pos_z > 0 and pos_z <= @exit", inplace=True)
            if df1.empty:
                    m = df
            else:
                # Remove df1 from df
                # m = df.merge(df1, on=['part_number'])
                m = pd.concat([df, df1]).drop_duplicates(keep=False)

        else:
            m = df

        bdsim_beam = pd.DataFrame(columns=['X', 'PX', 'Y', 'PY', 'P', 'E', 'ParentID', 'PDG_ID', 'Weight'])
        bdsim_beam['X'] = m[element.name + ".x"].values
        bdsim_beam['PX'] = m[element.name + ".xp"].values
        bdsim_beam['Y'] = m[element.name + ".y"].values
        bdsim_beam['PY'] = m[element.name + ".yp"].values
        bdsim_beam['E'] = m[element.name + ".energy"].values
        bdsim_beam['ParentID'] = m[element.name + ".parentID"].values
        bdsim_beam['PDG_ID'] = m[element.name + ".partID"].values
        bdsim_beam['Weight'] = m[element.name + ".weight"].values
        bdsim_beam.query('X == X', inplace=True)  # Remove None entries
        bdsim_beam.query("PDG_ID == 2212", inplace=True)
        bdsim_beam['E'] = 1000*bdsim_beam['E'] - physics.PROTON_MASS
        bdsim_beam['P'] = physics.energy_to_momentum(bdsim_beam['E'])

        # Query particles inside the aperture
        # TODO make a kwargs if we want to have an aperture (if using_collimators == TRUE)
        bdsim_beam["R"] = np.sqrt(bdsim_beam['X']**2 + bdsim_beam['Y']**2)

        if kwargs.get("with_aperture", True):
            if element['APERTYPE'] == 'CIRCLE':
                aperture = float(element['APERTURE'])
                bdsim_beam.query("R < @aperture", inplace=True)

            if element['APERTYPE'] == 'RECTANGLE':

                if element['TYPE'] == 'SLITS':
                    context = kwargs.get('context', {})
                    aper1 = context.get(f"w{element.name}X", 0.1)
                    aper2 = context.get(f"w{element.name}Y", 0.1)

                else:  # TODO Other case than bend ?
                    aperture = element['APERTURE'].split(",")
                    aper1 = 0.5*float(aperture[0])
                    aper2 = 0.5*float(aperture[1])

                bdsim_beam.query("abs(X) < @aper1 and abs(Y) < @aper2", inplace=True)

        tmp = beam_bdsim.BeamBdsim(bdsim_beam[['X', 'PX', 'Y', 'PY', 'P', 'E', 'ParentID', 'PDG_ID', 'Weight']])
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

    bd.track(bd_beam, p0, **kwargs)

    # Add options for Bdsim
    bd.set_options(options)

    # Run bdsim
    errors = bd.run(**kwargs).fatals

    # Create a new beamline to include the results
    l = line.line.copy()

    # TODO Make the extraction of the beam in MT : make with a kwargs
    # num_cores = multiprocessing.cpu_count()
    # results = Parallel(n_jobs=num_cores - 1)(delayed(extract_data)(i, line) for i, line in Data_simulation.iterrows())

    if kwargs.get("extract_beam", True):

        # Open the ROOT file and get the events
        f = ROOT.TFile(bd._outputname+'.root')
        evttree = f.Get("Event")
        # #
        # # # Add columns which contains datas
        print("WRITE BEAM FILE")
        l['BEAM'] = l.apply(lambda g: read_tracking(g, evttree, **kwargs), axis=1)
        #
    return beamline.Beamline(l)
