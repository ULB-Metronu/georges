import subprocess as sub
import re
import pandas as pd
import numpy as np
from ..simulator import Simulator
from ..simulator import SimulatorException
from ..lib.pybdsim import pybdsim
from .. import physics

INPUT_FILENAME = 'input.bdsim'

SUPPORTED_OPTIONS = [
            'beampipeRadius',
            'beampipeThickness',
            'beampipeMaterial',
            'circular',
            'elossHistoBinWidth',
            'eventNumberOffset',
            'hStyle',
            'vhRatio',
            'coilWidthFraction',
            'coilHeightFraction',
            'killNeutrinos',
            'ngenerate',
            'nturns',
            'outerDiameter',
            'physicsList',
            'printModuloFraction',
            'recreate',
            'recreateFileName',
            'startFromEvent',
            'seed',
            'seedStateFileName',
            'stopSecondaries',
            'useASCIISeedState',
            'writeSeedState',
            'aper1',
            'aper2',
            'aper3',
            'aper4',
            'dontSplitSBends',
            'ignoreLocalAperture',
            'checkOverlaps',
            'includeIronMagFields',
            'magnetGeometryType',
            'outerMaterial',
            'samplerDiameter',
            'sensitiveBeamlineComponents',
            'sensitiveBeamPipe',
            'vacuumMaterial',
            'vacuumPressure',
            'thinElementLength',
            'deltaChord',
            'deltaIntersection',
            'chordStepMinimum',
            'includeFringeFields',
            'integratorSet',
            'lengthSafety',
            'maximumEpsilonStep',
            'maximumStepLength',
            'maximumTrackingTime',
            'maximumTrackLength',
            'minimumEpsilonStep',
            'minimumRadiusOfCurvature',
            'deltaOneStep',
            'defaultBiasVacuum',
            'defaultBiasMaterial',
            'synchRadOn',
            'turnOnCerenkov',
            'defaultRangeCut',
            'prodCutPhotons',
            'prodCutElectrons',
            'prodCutPositrons',
            'prodCutProtons',
            'minimumKineticEnergy',
            'minimumRange',
            'storeTrajectories',
            'storeTrajectoryDepth',
            'storeTrajectoryEnergyThreshold',
            'storeTrajectoryParticle',
            'trajCutGTZ',
            'trajCutLTR',
            'nperfile',
            'nSegmentsPerCircle',
        ]


SUPPORTED_ELEMENTS = (
    'DRIFT',
    'GAP',
    'SBEND',
    'QUADRUPOLE',
    'SEXTUPOLE',
    'OCTUPOLE',
    'SLITS',
    'COLLIMATOR',
    'MARKER',
    'SOLIDS'
)


def sequence_to_bdsim(sequence, **kwargs):
    """Convert a pandas.DataFrame sequence onto a `pybdsim` input machine."""
    # Verify that all element types are supported
    if not sequence.apply(lambda e: e['TYPE'] in SUPPORTED_ELEMENTS, axis=1).all():
        raise BdsimException("Unsupported element type present in the sequence.")
    m = pybdsim.Builder.Machine()
    context = kwargs.get('context', {})
    for index, element in sequence.iterrows():
        if element['TYPE'] == 'DRIFT' and pd.notna(element['PIPE']):
            if pd.isna(element['APERTYPE']):
                m.AddDrift(index,
                           element['LENGTH']
                           )
            else:
                m.AddDrift(index,
                           element['LENGTH'],
                           apertureType=element['APERTYPE'],
                           aper1=(element['APERTURE'], 'm')
                           )
        if (element['TYPE'] == 'DRIFT' and element['PIPE'] is False) or element['TYPE'] == 'GAP':
            m.AddGap(index, element['LENGTH'])
        if element['TYPE'] == 'QUADRUPOLE':
            m.AddQuadrupole(index, element['LENGTH'], k1=context.get(f"{index}_K1", 0.0))
        if element['TYPE'] == 'SEXTUPOLE':
            m.AddSextupole(index, element['LENGTH'], k2=context.get(f"{index}_K2", 0.0))
        if element['TYPE'] == 'OCTUPOLE':
            m.AddSextupole(index, element['LENGTH'], k3=context.get(f"{index}_K3", 0.0))
        if element['TYPE'] == 'COLLIMATOR':
            m.AddECol(index,
                      element['LENGTH'],
                      xsize=float(element['APERTURE']),
                      ysize=float(element['APERTURE']),
                      material='Copper',
                      )
        if element['TYPE'] == "SLITS":
            m.AddRCol(index,
                      element['LENGTH'],
                      xsize=0.1,
                      ysize=0.1,
                      material="copper"
                     )
        if element['TYPE'] == "MARKER":
            m.AddMarker(index)
        if element['TYPE'] == "SOLIDS":
            m.AddGap(f"{index}_GAP", element['LENGTH'])
            m.AddPlacement(
                index,
                geometryFile=element['SOLIDS_FILE'],
                x=0.0,
                y=-0.229,
                z=element['AT_CENTER']-50/1000,
                phi=0.0,
                psi=np.deg2rad(30),
                theta=np.deg2rad(90)
            )
        if element['TYPE'] == "SBEND":
            m.AddDipole(
                index,
                'sbend',
                element['LENGTH'],
                angle=element['ANGLE'],
                e1=element['E1'],
                e2=element['E2']
            )
    m.AddSampler("all")
    return m


class BdsimException(Exception):
    """Exception raised for errors in the Twiss module."""

    def __init__(self, m):
        self.message = m


class BDSim(Simulator):
    """A Python wrapper around the BDSim software.

    Uses `pybdsim` .
    """

    EXECUTABLE_NAME = 'bdsim'

    def __init__(self, **kwargs):
        self._bdsim_machine = None
        self._bdsim_options = pybdsim.Options.Options()
        self._exec = BDSim.EXECUTABLE_NAME
        super().__init__(**kwargs)

    def _attach(self, beamline):
        super()._attach(beamline)
        if beamline.length is None or pd.isnull(beamline.length):
            raise SimulatorException("Beamline length not defined.")
        self._bdsim_machine = sequence_to_bdsim(beamline.line)

    def run(self, **kwargs):
        """Run bdsim as a subprocess."""

        if self._get_exec() is None:
            raise BdsimException("Can't run BDSim if no valid path and executable are defined.")

        p = sub.Popen(f"{self._get_exec()} --file={INPUT_FILENAME}",
                      stdin=sub.PIPE,
                      stdout=sub.PIPE,
                      stderr=sub.STDOUT,
                      cwd=".",
                      shell=True
                      )
        self._output = p.communicate()[0].decode()
        self._warnings = [line for line in self._output.split('\n') if re.search('warning|fatal', line)]
        self._fatals = [line for line in self._output.split('\n') if re.search('fatal', line)]
        self._last_context = kwargs.get("context", {})
        if kwargs.get('debug', False):
            print(self._output)
        return self

    def add_options(self, **kwargs):
        for k, v in kwargs.items():
            self.add_option(k, v)

    def add_option(self, key, value):
        if key not in SUPPORTED_OPTIONS:
            raise BdsimException("Invalid option.")
        # Ugly way to get the first letter capitalized...
        getattr(self._bdsim_options, f"Set{key[0].upper() + key[1:] }")(value)

    def beam_from_file(self, beam_file, **kwargs):
        b = pybdsim.Beam.Beam(
            particletype=kwargs.get('particletype', 'proton'),
            energy=kwargs.get('energy', physics.PROTON_MASS/1000),
            distrtype='userfile',
            distrFile=f"\"{beam_file}\"",
            distrFileFormat='"x[m]:xp[rad]:y[m]:yp[rad]:E[MeV]"'
        )
        self._bdsim_machine.AddBeam(b)
