:tocdepth: 1

PTW
===


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
--------------------

.. jupyter-execute::

    Bragg_peak_data = pd.read_csv(os.path.join(path,'Bragg_peak_library.csv'))
    z_values = np.arange(0,100.25,0.25)
    z_values = (z_values[0:-1] + (z_values[1:] - z_values[0:-1])[0]/2)[0:155]
    Bragg_peak_data['z'] = z_values

Compute the SOBP
----------------

.. jupyter-execute::

    sobp_analysis = SpreadOutBraggPeakAnalysis(dose_data=Bragg_peak_data.drop(columns='z').T.iloc[0:],
                                            method='scipy.optimize',
                                            z_axis=z_values[0:155])

.. jupyter-execute::

    sobp_analysis.compute_weights()

.. jupyter-execute::

    sobp_analysis.view_sobp(with_pristine_peaks=True)

Load lateral profile
--------------------

.. jupyter-execute::

    profile_df = pd.read_csv(os.path.join(path,'lateral_profile.csv'))
    lp_analysis = LateralProfileAnalysis(dose_profile=profile_df['dose'],
    positions=profile_df['x'])

.. jupyter-execute::

    lp_analysis.get_field_size()

.. jupyter-execute::

    lp_analysis.get_penumbra()

Module structure
----------------

.. automodapi:: georges.ptw
    :no-heading: