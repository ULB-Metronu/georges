import numpy as np
from georges import Kinematics
from georges import ureg as _
from georges import PlacementSequence, Element
from georges.manzoni import Input, Beam, SigmaObserver, BeamObserver, track
from georges.madx import MadX

# Create a test sequence by manually placing all elements
sequence = PlacementSequence()
sequence.place(Element.Quadrupole('Q1', L=1 * _.m, K1=1.0 / _.m**2), at=1.0 * _.m)
sequence.place(Element.Quadrupole('Q2', L=1 * _.m, K1=-1.0 / _.m**2), at=3.0 * _.m)
sequence.place(Element.Quadrupole('Q3', L=1 * _.m, K1=1.0 / _.m**2), at=5.0 * _.m)
sequence.place(Element.Marker('M1'), at=7.0 * _.m)
sequence.expand()  # Expand the sequence by creating implicit drifts

# Convert the sequence onto a Manzoni input
manzoni_input = Input.from_sequence(sequence)

# Create a Manzoni beam defined with a multi-Gaussian distribution
covariance = np.eye(6) / 200.0**2
covariance[4, 4] = 0.0
covariance[5, 5] = 0.0
distribution = np.random.multivariate_normal(np.zeros(6), covariance, int(1e3))
beam = Beam(kinematics=Kinematics(0.5), distribution=distribution)

# Create an observer to obtain the beam sizes
observer = SigmaObserver()

# The tracking is done here
track(beamline=manzoni_input, beam=beam, observer=observer)

# Convert the results to a nice DataFrame
results = observer.to_df()

# Test a few things
assert np.isclose(results.iloc[0]['BEAM_IN_X'], distribution[:, 0].std())
assert np.isclose(results.iloc[0]['BEAM_IN_XP'], distribution[:, 1].std())
assert np.isclose(results.iloc[0]['BEAM_IN_Y'], distribution[:, 2].std())
assert np.isclose(results.iloc[0]['BEAM_IN_YP'], distribution[:, 3].std())

# Now let's compare with MAD-X itself
# This is really cool: the sequence is automatically converted and sent to MAD-X
m = MadX(sequence=sequence, kinematics=Kinematics(0.5), stdout=False, command_log=None)
m.command.track(DELTAP=0.0, ONEPASS=True, ONETABLE=True, QUANTUM=False)
for i in range(distribution.shape[0]):
    m.command.start(x=distribution[i, 0],
                    px=distribution[i, 1],
                    y=distribution[i, 2],
                    py=distribution[i, 3],
                    t=0.0,
                    pt=distribution[i, 4])
m.command.run(turns=1)
m.command.endtrack()

results_madx = m.table.trackone.dframe()

beam_observer = BeamObserver()
track(beamline=manzoni_input, beam=beam, observer=beam_observer)
results_manzoni = beam_observer.to_df()

np.all(np.isclose(
    results_madx.query("turn == 1.0")[['x', 'px', 'y', 'py']].values,
    results_manzoni.iloc[-1]['BEAM_OUT'][:, 0:4]
))
