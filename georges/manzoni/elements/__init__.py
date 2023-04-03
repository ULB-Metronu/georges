"""
The `elements` submodule contains the physical implementation of the different beamline elements,
namely:

- Drifts
- Dipoles
- Quadrupoles
- Sextupoles
- Multipoles (magnets with a combined dipolar, quadrupolar and sextupolar function),
- Collimators (rectangular, circular and elliptical);
- Scatterers
    thin materials that interact with the particles and only modify the transverse angle along the two transverse axes
- Degraders
    thick materials which combine particle interaction with a simple drift type propagation along the  material
"""
from .collimators import CircularCollimator, Collimator, Dump, EllipticalCollimator, RectangularCollimator
from .elements import ManzoniElement
from .magnets import (
    Bend,
    DipEdge,
    Drift,
    Fringein,
    Fringeout,
    Gap,
    HKicker,
    Kicker,
    Marker,
    Matrix,
    Multipole,
    Quadrupole,
    RBend,
    SBend,
    Sextupole,
    SRotation,
    TKicker,
    VKicker,
)
from .scatterers import BeamStop, Degrader, Scatterer
