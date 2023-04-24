import os
import pandas as pd
from numpy.testing import assert_almost_equal

import georges

from georges.ptw import LateralProfileAnalysis


def test_lateral_profile():
    path = f"{os.path.join(georges.__path__[0], '../docs/ptw_files')}"
    profile_df = pd.read_csv(os.path.join(path, "lateral_profile.csv"))
    lp_analysis = LateralProfileAnalysis(dose_profile=profile_df["dose"], positions=profile_df["x"])

    assert_almost_equal(lp_analysis.get_field_size(), 29.327797, decimal=4)
    assert_almost_equal(lp_analysis.get_penumbra(), 0.7944795917355796, decimal=4)
    assert_almost_equal(lp_analysis.get_position_left(12), -14.949340467988621, decimal=4)
    assert_almost_equal(lp_analysis.get_position_right(1), 18.105823113508254, decimal=4)
    assert_almost_equal(lp_analysis.get_ur_size(), 26.14987903551072, decimal=4)
    assert_almost_equal(lp_analysis.get_ur_flatness(), 3.06115, decimal=4)
