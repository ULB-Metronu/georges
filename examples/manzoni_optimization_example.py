import numpy as np
import matplotlib.pyplot as plt
import dlib
from georges import Kinematics
from georges import ureg as _
from georges import PlacementSequence, Element
from georges.manzoni import Input, Beam, SigmaObserver, track

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
covariance = np.eye(5) / 200.0**2
distribution = np.random.multivariate_normal(np.zeros(5), covariance, int(1e5))
beam = Beam(kinematics=Kinematics(0.5), distribution=distribution)

# Create an observer to obtain the beam sizes
observer = SigmaObserver()


# Create the optimization routine
def parametric_tracking(kq1, kq2):
    manzoni_input.sequence[1].K1 = kq1 / _.m**2
    manzoni_input.sequence[3].K1 = kq2 / _.m**2
    manzoni_input.sequence[5].K1 = kq2 / _.m**2
    track(beamline=manzoni_input, beam=beam, observer=observer)
    if np.isnan(observer.data[-1][1]) or np.isnan(observer.data[-1][3]):
        return 1000.0
    return observer.data[-1][1]**2 + observer.data[-1][3]


iters = []
v = []
for i in [5, 10, 25, 50, 75, 100, 125]:
    x, y = dlib.find_min_global(parametric_tracking,
                                [0, -3.5],  # Lower bound constraints
                                [3.5, 0.0],  # Upper bound constraints
                                i)
    v.append(y)
    print(y)
    iters.append(i)

plt.plot(iters, v, '*-')
plt.show()
