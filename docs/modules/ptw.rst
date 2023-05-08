:tocdepth: 2

***
PTW
***

The PTW module includes several classes that allow for direct assessment of various clinical parameters in charged particle therapy. These parameters include distal fall-off, range, transverse and in-depth flatness, and lateral penumbra. The following classes are implemented:

* **BraggPeakAnalysis**: This class processes a normalized 1D depth dose profile, or a pristine Bragg peak, to extract key information such as the maximum dose position, the distal range at a given percentage of the maximum, and the distal fall-off (DFO). The user must provide both the dose data and the corresponding positions in depth.

* **SpreadOutBraggPeakAnalysis**: This class takes a set of Bragg peaks, or a Bragg peak library, as input and computes the relative weights required to obtain a uniform depth dose profile, known as the SOBP. Additionally, this class processes the resulting SOBP and provides the flatness, DFO, and distal range at a given percentage of the maximum to the user.

* **LateralProfileAnalysis**: This class takes a normalized 1D transverse dose profile as input and computes various parameters such as the field size, the uniform region (defined as 80% of the field size), the transverse flatness of the uniform region, and the lateral penumbra at the left and right of this region. The user must provide both the dose data and the corresponding transverse positions.

The module also includes two classes that help optimize spot spacing in pencil beam scanning treatment:

* **RegularSpotScanning**: This class provides a method to calculate the required spot spacing for a regular grid irradiation scheme in order to obtain a two-dimensional uniform dose deposition profile for a given spot width (1 sigma) and a targeted circular field. The user inputs the half value of the field size, the desired number of spots along each axis of the field, and the standard deviation (1 sigma) of the beam. The required spot spacing to achieve a 2D dose uniformity of at least 98% is directly outputted.

* **ContourSpotScanning**: This class works similarly to RegularSpotScanning but uses a circular, contour-based irradiation scheme with a central spot placed at the center of the field. The user can choose whether to impose the irradiation radius. The output of the calculation includes the relative weight of the contour spots compared to the central spot and the angle spacing between these spots.

Example
#######


.. jupyter-execute::
    :hide-output:

    import os
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt
    import georges

    from georges.ptw import SpreadOutBraggPeakAnalysis, \
                            BraggPeakAnalysis, \
                            LateralProfileAnalysis


.. jupyter-execute::
    :hide-output:
    :hide-code:

    path = f"{os.path.join(georges.__path__[0], '../docs/ptw_files')}"

Load Bragg peak data
""""""""""""""""""""

.. jupyter-execute::

    Bragg_peak_data = pd.read_csv(os.path.join(path,'Bragg_peak_library.csv'))
    z_values = np.arange(0,100.25,0.25)
    z_values = (z_values[0:-1] + (z_values[1:] - z_values[0:-1])[0]/2)[0:155]
    Bragg_peak_data['z'] = z_values

Compute the SOBP
""""""""""""""""

.. jupyter-execute::

    sobp_analysis = SpreadOutBraggPeakAnalysis(dose_data=Bragg_peak_data.drop(columns='z').T.iloc[0:],
                                            method='scipy.optimize',
                                            z_axis=z_values[0:155])

.. jupyter-execute::

    sobp_analysis.compute_weights()

.. jupyter-execute::

    sobp_analysis.view_sobp(with_pristine_peaks=True)

Load lateral profile
""""""""""""""""""""

.. jupyter-execute::

    profile_df = pd.read_csv(os.path.join(path,'lateral_profile.csv'))
    lp_analysis = LateralProfileAnalysis(dose_profile=profile_df["dose"], positions=profile_df["x"])

.. jupyter-execute::

    lp_analysis.get_field_size()

.. jupyter-execute::

    lp_analysis.get_penumbra()

API
###

.. automodapi:: georges.ptw
    :no-heading: