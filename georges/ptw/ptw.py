"""
TODO
"""
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as _np
import pandas as pd
import scipy
from numpy.polynomial import Polynomial
from scipy.interpolate import interp1d
from scipy.optimize import Bounds, minimize, newton_krylov


def compute_dvh(dose_data, voxel_volume):
    """
    Args:
        dose_data: 3D-array of the dose values.
        voxel_volume: The volume of a unique voxel. We consider the same size for all voxels.

    Returns:
        dvh_dataframe: The dvh stored into a 2 columns dataframe

    """
    dvh_histogram = plt.hist(
        dose_data.flatten(),
        bins=100,
        cumulative=-1,
        density=True,
    )

    dvh_dataframe = pd.DataFrame(
        columns=[
            "volume",
            "dose_value",
        ],
    )

    dvh_dataframe["volume"] = dvh_histogram[0]
    dvh_dataframe["dose_value"] = voxel_volume * (
        dvh_histogram[1][:-1] + (dvh_histogram[1][1] - dvh_histogram[1][0]) / 2
    )

    return dvh_dataframe


class BraggPeakException(Exception):
    def __init__(self, m: str = "") -> None:  # pragma: no cover
        self.message = m


class BraggPeakAnalysis:
    """
    This class processes a normalized 1D depth dose profile, or a pristine Bragg peak, to extract key information such
    as the maximum dose position, the distal range at a given percentage of the maximum, and the distal fall-off (DFO).
    The user must provide both the dose data and the corresponding positions in depth.
    """

    def __init__(self, bp: pd.DataFrame, method: str, low_dose=1, high_dose=70):
        self.data = bp
        self.low_bound = low_dose
        self.high_bound = high_dose
        self.coefficients = None
        self.bp_interval = None
        self.method = method

    def get_distal_interval(self):
        idx = _np.where(_np.isclose(self.data[self.data.columns[1]], 1e2))
        distal_interval = self.data[idx[0][0] :]

        return distal_interval

    def get_bp_interval(self):
        idx1 = _np.where(self.data[self.data.columns[1]] > self.high_bound)
        idx2 = _np.where(self.data[self.data.columns[1]] < self.low_bound)

        idx_min = idx1[0][0]
        idx_max = idx2[0][0]
        bp_interval = self.data[idx_min:idx_max]

        return bp_interval

    def compute_percentage(self, x):  # Gives dose percentage in the distal region for a given z position
        f = interp1d(
            self.get_distal_interval()[self.get_distal_interval().columns[0]].values,
            self.get_distal_interval()[self.get_distal_interval().columns[1]].values,
            kind=2,
            bounds_error=False,
        )

        return f(x)

    def fitting_function(self, x):
        coeff = self.coefficients.convert().coef
        return (
            coeff[0]
            + coeff[1] * x
            + coeff[2] * x**2
            + coeff[3] * x**3
            + coeff[4] * x**4
            + coeff[5] * x**5
            + coeff[6] * x**6
            + coeff[7] * x**7
            + coeff[8] * x**8
        )

    def fit_bp(self):
        """
        This function fits a Bragg peak using a 8th order polynom. Based on reference xxx.

        Args:
        Returns:

        """

        bp_interval = self.get_bp_interval()
        bp_interval.reset_index(drop=True, inplace=True)

        x = _np.float32(bp_interval[bp_interval.columns[0]])
        y = _np.float32(bp_interval[bp_interval.columns[1]])

        self.coefficients = Polynomial.fit(
            x=x,
            y=y,
            deg=8,
            full=False,
        )

        self.bp_interval = bp_interval[bp_interval.columns[0]].to_numpy()

        return self.bp_interval, self.coefficients, self.fitting_function(self.bp_interval)

    def get_coefficients(self) -> _np.array:
        return self.coefficients.convert().coef

    def get_xrange(self, val):
        if self.method == "polynomial_fit":
            rval = (self.fit_bp()[1] - val).roots()
            real_val = rval[rval.imag == 0].real
            var = real_val[
                _np.where(
                    _np.logical_and(
                        real_val >= self.get_distal_interval()[self.get_distal_interval().columns[0]].values[0],
                        real_val <= self.get_distal_interval()[self.get_distal_interval().columns[0]].values[-1],
                    ),
                )
            ]
            # check if we are in the distal interval of the BP
            if var.size:
                return _np.max(var)
            return 0

        elif self.method == "scipy.optimize":
            nearest_val = self.get_distal_interval()[self.get_distal_interval().columns[0]][
                _np.abs(self.get_distal_interval()[self.get_distal_interval().columns[1]] - val).idxmin()
            ]

            var = newton_krylov(
                lambda x: _np.abs(self.compute_percentage(x=x) - val),
                xin=nearest_val,
                maxiter=10000,
                f_tol=1e-12,
                x_tol=1e-13,
            )

            return var

        else:
            raise BraggPeakException(
                "Please choose a valid computation method (must be either 'polynomial_fit' or 'scipy.optimize')",
            )

    def get_r90(self) -> float:
        return self.get_xrange(90)

    def get_r80(self) -> float:
        return self.get_xrange(80)

    def get_r20(self) -> float:
        return self.get_xrange(20)

    def get_r10(self) -> float:
        return self.get_xrange(10)

    def get_maximum(self) -> float:
        roots = self.fit_bp()[1].deriv().roots()
        real_roots = roots[roots.imag == 0].real
        var = real_roots[
            _np.where(
                _np.logical_and(
                    real_roots >= self.get_distal_interval()[self.get_distal_interval().columns[0]].values[0],
                    real_roots <= self.get_distal_interval()[self.get_distal_interval().columns[0]].values[-1],
                ),
            )
        ]
        fit_var = self.fitting_function(var)
        var = var[_np.where(fit_var == _np.max(fit_var))]  # take the value where the fit is maximum (global maximum)
        return var[0]

    def compute_distal_fall_off(self) -> List[float]:
        return [self.get_r10() - self.get_r90(), self.get_r20() - self.get_r80()]

    def compute_dfo_90_10(self) -> float:
        return self.compute_distal_fall_off()[0]

    def compute_dfo_80_20(self) -> float:
        return self.compute_distal_fall_off()[1]


class SpreadOutBraggPeakAnalysis:
    """
    This class takes a set of Bragg peaks, or a Bragg peak library, as input and computes the relative weights required
    to obtain a uniform depth dose profile, known as the SOBP. Additionally, this class processes the resulting SOBP
    and provides the flatness, DFO, and distal range at a given percentage of the maximum to the user.
    """

    def __init__(
        self,
        dose_data: Optional[pd.DataFrame] = None,
        method: str = scipy.optimize,
        z_axis: Optional[_np.array] = None,
        color: str = "blue",
        str_on_legend: str = "",
        adjust_last_peak: float = 1.0,
    ):
        self.dose_data = dose_data
        self.method = method
        self.z_axis = z_axis
        self.modul_type = "Full"
        self.color = color
        self.str_on_legend = str_on_legend
        self.adjust_last_peak = adjust_last_peak

    def get_library_max_ranges(self):
        max_ranges = _np.zeros(self.dose_data.shape[0])

        for i in range(self.dose_data.shape[0]):
            normalized_bragg_peak = 1e2 * self.dose_data.iloc[i, :] / _np.max(self.dose_data.iloc[i, :])
            bp_analysis = BraggPeakAnalysis(
                bp=pd.DataFrame(
                    {
                        "centers": self.z_axis,
                        "dose": normalized_bragg_peak,
                    },
                ),
                method="scipy.optimize",
            )

            max_ranges[i] = bp_analysis.get_xrange(100)

        return max_ranges

    def sobp_data(self) -> _np.array:
        sobp_data = _np.zeros(
            (
                self.dose_data.shape[0],
                self.dose_data.shape[0],
            ),
        )  # Array to store maximum dose values of each Bragg Peak

        max_ranges = self.get_library_max_ranges()
        for i in range(self.dose_data.shape[0]):
            for j in range(self.dose_data.shape[0]):
                bp_analysis = BraggPeakAnalysis(
                    bp=pd.DataFrame(
                        {
                            "z": self.z_axis,
                            "dose": self.dose_data.iloc[j, :],
                        },
                    ),
                    method="scipy.optimize",
                )
                f = interp1d(
                    bp_analysis.data[bp_analysis.data.columns[0]].values,
                    bp_analysis.data[bp_analysis.data.columns[1]].values,
                    kind=2,
                    bounds_error=False,
                )
                sobp_data[i, j] = f(max_ranges[i])

        return sobp_data

    def compute_weights(self) -> _np.array:
        a_matrix = self.sobp_data() / self.sobp_data().max()  # Normalize all data to max 1
        goal_dose_values = _np.full(a_matrix.shape[0], 1)

        if self.method == "np.linalg.solve":  # This method may give negative weights !
            return _np.linalg.solve(a_matrix, goal_dose_values)

        elif self.method == "scipy.optimize":
            b = goal_dose_values
            n = len(b)
            fun = lambda x: _np.linalg.norm(_np.dot(a_matrix, x) - b)  # noqa: E731
            sol = minimize(
                fun,
                _np.zeros(n),
                method="L-BFGS-B",
                bounds=[(0.0, None) for _ in range(n)],
            )

            return sol["x"]

    def get_final_weights(self):
        if self.adjust_last_peak is not None:
            weights = self.compute_weights()
            weights[0] *= self.adjust_last_peak

        return weights

    def compute_sobp_profile(self) -> _np.array:
        weights = self.get_final_weights()

        weights = weights.reshape(1, self.dose_data.shape[0])
        weighted_dose_data = _np.matmul(
            weights,
            self.dose_data.values,
        )

        sobp_profile = weighted_dose_data.reshape(self.dose_data.shape[1])

        return 1e2 * sobp_profile / sobp_profile.max()

    def view_sobp(self, with_pristine_peaks=False):
        def plot(x, y):
            plt.plot(
                x,
                y,
                linestyle="dashed",
                linewidth=0.7,
                color="r",
                marker="*",
                markersize=4,
            )

        if with_pristine_peaks:
            weighted_dose = self.dose_data.multiply(self.compute_weights(), axis=0)
            weighted_dose.apply(
                lambda e: plot(
                    x=self.z_axis,
                    y=1e2 * e.values / self.sobp_data().max(),
                ),
                axis="columns",
            )

        plt.plot(
            self.z_axis,
            self.compute_sobp_profile(),
            linestyle="dashed",
            linewidth=0.8,
            color=self.color,
            marker="*",
            markersize=5,
            label=self.str_on_legend,
        )

        plt.xlabel("Depth in water (mm)")
        plt.ylabel("Normalized dose (\\%)")
        plt.legend(loc=(0.15, 1.01), fontsize=20, ncol=2, frameon=False)
        plt.xlim(0, _np.max(self.z_axis))

    def recompute_sobp_for_another_range(self, bp_to_leave):
        if bp_to_leave == 0:
            new_weighted_dose = _np.matmul(
                self.compute_weights().reshape(1, self.dose_data.shape[0]),
                self.dose_data.values[0:, :],
            )
        else:
            new_shape = self.dose_data.shape[0] - bp_to_leave
            weights = self.compute_weights()[0:new_shape]

            new_weighted_dose = _np.matmul(
                weights.reshape(1, new_shape),
                self.dose_data.values[bp_to_leave:, :],
            )

        sobp_data = new_weighted_dose.reshape(self.dose_data.shape[1])
        sobp_data = 1e2 * sobp_data / sobp_data.max()

        return sobp_data

    def recompute_sobp_for_another_width(self, bp_to_leave):
        new_shape = self.dose_data.shape[0] - bp_to_leave
        weights = self.compute_weights()[0:new_shape]
        new_weighted_dose = _np.matmul(
            weights.reshape(1, new_shape),
            self.dose_data.values[0:new_shape, :],
        )

        sobp_data = new_weighted_dose.reshape(self.dose_data.shape[1])
        sobp_data = 1e2 * sobp_data / sobp_data.max()

        return sobp_data

    def compute_ranges(self):
        sobp_array = self.compute_sobp_profile()
        idxmax = _np.where(sobp_array < 1)[0][0]  # Select the dose values higher than 1%

        idx = 0
        for i in range(_np.flip(sobp_array[0:idxmax]).shape[0]):
            if _np.flip(sobp_array[0:idxmax])[i + 1] > _np.flip(sobp_array[0:idxmax])[i]:
                idx += 1
            else:
                break

        flipped_data = _np.flip(sobp_array[idxmax - idx : idxmax])  # Select the ramp

        def compute_range(data, r, z_axis):
            idx_xin = pd.DataFrame(_np.abs(data - r)).idxmin()
            nearest_val = _np.flip(z_axis[idxmax - idx : idxmax])[idx_xin[0]]
            fct = interp1d(
                _np.flip(z_axis[idxmax - idx : idxmax]),
                data,
                kind=2,
                bounds_error=False,
            )

            res = minimize(
                fun=lambda x: _np.abs(fct(x=x) - r),
                bounds=Bounds(
                    lb=nearest_val - 0.5,
                    ub=nearest_val + 0.5,
                ),
                x0=nearest_val,
                method="trust-constr",
                options={
                    "xtol": 1e-8,
                    "verbose": 0,
                    "maxiter": 1e5,
                },
            )
            return res

        r_10 = compute_range(flipped_data, 10, self.z_axis).x[0]
        r_20 = compute_range(flipped_data, 20, self.z_axis).x[0]
        r_80 = compute_range(flipped_data, 80, self.z_axis).x[0]
        r_90 = compute_range(flipped_data, 90, self.z_axis).x[0]

        return r_10, r_20, r_80, r_90

    def get_sobp_r10(self):
        return self.compute_ranges()[0]

    def get_sobp_r20(self):
        return self.compute_ranges()[1]

    def get_sobp_r80(self):
        return self.compute_ranges()[2]

    def get_sobp_r90(self):
        return self.compute_ranges()[3]

    def get_sobp_flatness(self):
        sobp_array = self.compute_sobp_profile()
        idx_98 = _np.where(_np.isclose(sobp_array, 98, atol=1.5))[0]
        max_dose = _np.max(sobp_array[idx_98[0] : idx_98[-1]])
        min_dose = _np.min(sobp_array[idx_98[0] : idx_98[-1]])

        flatness = 1e2 * (max_dose - min_dose) / (max_dose + min_dose)

        return flatness

    def get_ur_min_dose(self):
        sobp_array = self.compute_sobp_profile()
        idx_98 = _np.where(_np.isclose(sobp_array, 98, atol=1.5))[0]
        min_dose = _np.min(sobp_array[idx_98[0] : idx_98[-1]])

        return min_dose

    def get_ur_max_dose(self):
        sobp_array = self.compute_sobp_profile()
        idx_98 = _np.where(_np.isclose(sobp_array, 98, atol=1.5))[0]
        max_dose = _np.max(sobp_array[idx_98[0] : idx_98[-1]])

        return max_dose


class LateralProfileAnalysis:
    """
    This class takes a normalized 1D transverse dose profile as input and computes various parameters such as the field
    size, the uniform region (defined as 80% of the field size), the transverse flatness of the uniform region, and the
    lateral penumbra at the left and right of this region. The user must provide both the dose data and the
    corresponding transverse positions.
    """

    def __init__(self, dose_profile: _np.array, positions: _np.array):
        self.dose_profile = dose_profile
        self.positions = positions

    def set_data(self):
        idxmax = _np.where(self.positions == self.positions[_np.int(self.positions.shape[0] / 2)])[0][0]
        positions_left = self.positions[0:idxmax]
        positions_right = self.positions[idxmax:]
        dose_left = self.dose_profile[0:idxmax]
        dose_right = self.dose_profile[idxmax:]

        return idxmax, positions_left, dose_left, positions_right, dose_right

    def define_f(self):
        f_left = interp1d(
            self.set_data()[1],
            self.set_data()[2],
            kind=1,
            bounds_error=False,
        )
        f_right = interp1d(
            self.set_data()[3],
            self.set_data()[4],
            kind=1,
            bounds_error=False,
        )
        return f_left, f_right

    def get_position_left(self, percentage):
        nearest_val = self.set_data()[1][pd.DataFrame(_np.abs(self.set_data()[2] - percentage)).idxmin()]
        position = minimize(
            fun=lambda x: _np.abs(self.define_f()[0](x=x) - percentage),
            bounds=Bounds(
                lb=nearest_val - 3,
                ub=nearest_val + 3,
            ),
            x0=nearest_val,
            method="trust-constr",
            options={
                "xtol": 1e-8,
                "verbose": 0,
                "maxiter": 1e5,
            },
        )
        return position.x[0]

    def get_position_right(self, percentage):
        nearest_val = self.set_data()[3][pd.DataFrame(_np.abs(self.set_data()[4] - percentage)).idxmin()]
        position = minimize(
            fun=lambda x: _np.abs(self.define_f()[1](x=x) - percentage),
            bounds=Bounds(
                lb=nearest_val - 3,
                ub=nearest_val + 3,
            ),
            x0=nearest_val,
            method="trust-constr",
            options={
                "xtol": 1e-8,
                "verbose": 0,
                "maxiter": 1e5,
            },
        )
        return position.x[0]

    def get_p_20_left(self):
        return self.get_position_left(20)

    def get_p_50_left(self):
        return self.get_position_left(50)

    def get_p_80_left(self):
        return self.get_position_left(80)

    def get_p_20_right(self):
        return self.get_position_right(20)

    def get_p_50_right(self):
        return self.get_position_right(50)

    def get_p_80_right(self):
        return self.get_position_right(80)

    def get_field_size(self):
        return self.get_p_50_right() - self.get_p_50_left()

    def get_ur_left(self):
        return self.get_p_50_left() + 2 * self.get_penumbra_left()

    def get_ur_right(self):
        return self.get_p_50_right() - 2 * self.get_penumbra_right()

    def get_ur_size(self):
        return self.get_ur_right() - self.get_ur_left()

    def get_ur_min_dose(self):
        ur_idx_min = _np.where(self.positions > self.get_ur_left())[0][0]
        ur_idx_max = _np.where(self.positions > self.get_ur_right())[0][0]

        return _np.min(self.dose_profile[ur_idx_min:ur_idx_max])

    def get_ur_max_dose(self):
        ur_idx_min = _np.where(self.positions > self.get_ur_left())[0][0]
        ur_idx_max = _np.where(self.positions > self.get_ur_right())[0][0]

        return _np.max(self.dose_profile[ur_idx_min:ur_idx_max])

    def get_ur_flatness(self):
        return (
            1e2 * (self.get_ur_max_dose() - self.get_ur_min_dose()) / (self.get_ur_max_dose() + self.get_ur_min_dose())
        )

    def get_penumbra_left(self):
        return self.get_p_80_left() - self.get_p_20_left()

    def get_penumbra_right(self):
        return self.get_p_20_right() - self.get_p_80_right()

    def get_penumbra(self):
        return _np.mean([self.get_penumbra_left(), self.get_penumbra_right()])

    def plot(self):
        plt.vlines(
            x=[
                self.get_p_20_left(),
                self.get_p_50_left(),
                self.get_p_80_left(),
                self.get_p_20_right(),
                self.get_p_50_right(),
                self.get_p_80_right(),
            ],
            ymin=0,
            ymax=100,
            color="k",
            linestyle="dashed",
        )

        plt.hlines(
            y=[100, self.get_ur_min_dose(), self.get_ur_max_dose()],
            xmin=self.get_p_20_left(),
            xmax=self.get_p_20_right(),
            color="k",
            linestyle="dashed",
        )

        plt.text(
            s=f"Flatness: {_np.round(self.get_ur_flatness(), 2)} \\%",
            x=-self.get_field_size() / 5,
            y=85,
            fontsize=15,
            color="b",
        )
        plt.text(
            s=f"P$_L$: {_np.round(self.get_penumbra_left(), 2)} mm",
            x=self.get_p_50_left() - 2 * self.get_penumbra_left(),
            y=30,
            fontsize=15,
            color="b",
        )
        plt.text(
            s=f"P$_R$: {_np.round(self.get_penumbra_right(), 2)} mm",
            x=self.get_p_50_right() + self.get_penumbra_right(),
            y=30,
            fontsize=15,
            color="b",
        )

        plt.plot(
            self.positions,
            self.dose_profile,
            marker="o",
            markersize=1.5,
            linestyle="dashed",
            color="b",
        )

        plt.xlabel("Transverse position (mm)")
        plt.ylabel("Normalized dose (\\%)")

        plt.grid()


class GammaAnalysis:
    pass


class RegularSpotScanning:
    """
    This class provides a method to calculate the required spot spacing for a regular grid irradiation scheme in order
    to obtain a two-dimensional uniform dose deposition profile for a given spot width (1 sigma) and a targeted
    circular field. The user inputs the half value of the field size, the desired number of spots along each axis of
    the field, and the standard deviation (1 sigma) of the beam. The required spot spacing to achieve a 2D dose
    uniformity of at least 98% is directly outputted.
    """

    def __init__(self, sigma, fieldsize, n_spots_per_axis):
        self.sigma = sigma
        self.fieldsize = fieldsize
        self.n_spots_per_axis = n_spots_per_axis

    def double_gaussian_function(self, x, y, a, mu_x, mu_y):
        return (
            a
            * _np.exp(
                -((x - mu_x) ** 2) / (2 * self.sigma**2) - (y - mu_y) ** 2 / (2 * self.sigma**2),
            )
            / (_np.sqrt(2 * _np.pi) * self.sigma)
        )

    def compute_2d_scanned_profile(self, spot_space):
        positions = _np.arange(-50, 51, 1)
        x, y = _np.meshgrid(positions, positions)
        shape = x.shape

        dose_df = pd.DataFrame(
            {
                "x": x.reshape(shape[0] * shape[1]),
                "y": y.reshape(shape[0] * shape[1]),
            },
        )

        for i in range(self.n_spots_per_axis):
            for j in range(self.n_spots_per_axis):
                dose_df[f"dose_{i}_{j}"] = self.double_gaussian_function(
                    x=dose_df["x"],
                    y=dose_df["y"],
                    a=100,
                    mu_x=(-(self.n_spots_per_axis - 1) / 2 + i) * spot_space,
                    mu_y=(-(self.n_spots_per_axis - 1) / 2 + j) * spot_space,
                )

        total_dose = dose_df.drop(columns=["x", "y"]).sum(axis=1)
        dose_df["total_dose"] = total_dose

        return dose_df

    def optimize_regular_grid_spots_placement(self):
        def compute_2d_profile_flatness(x):
            df = self.compute_2d_scanned_profile(spot_space=x)
            profile = df.query(f"x**2+y**2<={self.fieldsize**2}")["total_dose"]

            flatness = 1e2 * (_np.max(profile) - _np.min(profile)) / (_np.max(profile) + _np.min(profile))

            return flatness

        optim = minimize(
            fun=compute_2d_profile_flatness,
            bounds=Bounds(
                lb=0,
                ub=5 * self.sigma,
            ),
            x0=3 * self.sigma,
            method="trust-constr",
            options={
                "verbose": 1,
                "xtol": 1e-20,
                "maxiter": 1e6,
            },
        )

        return optim


class ContourSpotScanning:
    """
    This class works similarly to RegularSpotScanning but uses a circular, contour-based irradiation scheme with a
    central spot placed at the center of the field. The user can choose whether to impose the irradiation radius. The
    output of the calculation includes the relative weight of the contour spots compared to the central spot and the
    angle spacing between these spots.
    """

    def __init__(self, sigma, fieldsize, desired_angle, shoot_on_aperture: bool = True, angle_imposed: bool = False):
        self.sigma = sigma
        self.fieldsize = fieldsize
        self.desired_angle = desired_angle
        self.shoot_on_aperture = shoot_on_aperture
        self.angle_imposed = angle_imposed

    def double_gaussian_function(self, x, y, a, mu_x, mu_y):
        return (
            a
            * _np.exp(
                -((x - mu_x) ** 2) / (2 * self.sigma**2) - (y - mu_y) ** 2 / (2 * self.sigma**2),
            )
            / (_np.sqrt(2 * _np.pi) * self.sigma)
        )

    def compute_2d_contour_scanned_profile(self, angle, weight):
        positions = _np.arange(-50, 51, 1)
        x, y = _np.meshgrid(
            positions,
            positions,
        )

        n_contour_spots = _np.int(2 * _np.pi / angle)
        dose_df = pd.DataFrame(
            {
                "x": x.reshape(x.shape[0] * x.shape[1]),
                "y": y.reshape(x.shape[0] * x.shape[1]),
            },
        )

        dose_df["dose_central_spot"] = self.double_gaussian_function(
            x=dose_df["x"],
            y=dose_df["y"],
            a=100,
            mu_x=0.0,
            mu_y=0.0,
        )
        for i in range(n_contour_spots):
            tmp_theta = -_np.pi + i * angle
            dose_df[f"dose_{i}"] = (
                self.double_gaussian_function(
                    x=dose_df["x"],
                    y=dose_df["y"],
                    a=100,
                    mu_x=self.fieldsize * _np.cos(tmp_theta),
                    mu_y=self.fieldsize * _np.sin(tmp_theta),
                )
                * weight
            )

        total_dose = dose_df.drop(columns=["x", "y"]).sum(axis=1)
        dose_df["total_dose"] = total_dose

        return dose_df

    def compute_2d_contour_scanned_profile_angle_imposed(self, weight):
        positions = _np.arange(-50, 51, 1)
        x, y = _np.meshgrid(
            positions,
            positions,
        )

        n_contour_spots = _np.int(2 * _np.pi / self.desired_angle)
        dose_df = pd.DataFrame(
            {
                "x": x.reshape(x.shape[0] * x.shape[1]),
                "y": y.reshape(x.shape[0] * x.shape[1]),
            },
        )

        dose_df["dose_central_spot"] = self.double_gaussian_function(
            x=dose_df["x"],
            y=dose_df["y"],
            a=100,
            mu_x=0.0,
            mu_y=0.0,
        )
        for i in range(n_contour_spots):
            tmp_theta = -_np.pi + i * self.desired_angle
            dose_df[f"dose_{i}"] = (
                self.double_gaussian_function(
                    x=dose_df["x"],
                    y=dose_df["y"],
                    a=100,
                    mu_x=self.fieldsize * _np.cos(tmp_theta),
                    mu_y=self.fieldsize * _np.sin(tmp_theta),
                )
                * weight
            )

        total_dose = dose_df.drop(columns=["x", "y"]).sum(axis=1)
        dose_df["total_dose"] = total_dose

        return dose_df

    def compute_2d_contour_scanned_profile_radius_not_imposed(self, angle, weight, shoot_at):
        positions = _np.arange(-50, 51, 1)
        x, y = _np.meshgrid(
            positions,
            positions,
        )

        n_contour_spots = _np.int(2 * _np.pi / angle)
        dose_df = pd.DataFrame(
            {
                "x": x.reshape(x.shape[0] * x.shape[1]),
                "y": y.reshape(x.shape[0] * x.shape[1]),
            },
        )

        dose_df["dose_central_spot"] = self.double_gaussian_function(
            x=dose_df["x"],
            y=dose_df["y"],
            a=100,
            mu_x=0.0,
            mu_y=0.0,
        )
        for i in range(n_contour_spots):
            tmp_theta = -_np.pi + i * angle
            dose_df[f"dose_{i}"] = (
                self.double_gaussian_function(
                    x=dose_df["x"],
                    y=dose_df["y"],
                    a=100,
                    mu_x=shoot_at * _np.cos(tmp_theta),
                    mu_y=shoot_at * _np.sin(tmp_theta),
                )
                * weight
            )

        total_dose = dose_df.drop(columns=["x", "y"]).sum(axis=1)
        dose_df["total_dose"] = total_dose

        return dose_df

    def optimize_contour_spots(self):
        def compute_contour_profile_flatness(x):
            df = self.compute_2d_contour_scanned_profile(
                angle=x[0],
                weight=x[1],
            )
            profile = df.query("x**2+y**2<=@self.fieldsize**2")["total_dose"]
            flatness = 1e2 * (_np.max(profile) - _np.min(profile)) / (_np.max(profile) + _np.min(profile))

            return flatness

        def compute_contour_profile_flatness_angle_imposed(x):
            df = self.compute_2d_contour_scanned_profile_angle_imposed(weight=x)
            profile = df.query("x**2+y**2<=@self.fieldsize**2")["total_dose"]
            flatness = 1e2 * (_np.max(profile) - _np.min(profile)) / (_np.max(profile) + _np.min(profile))

            return flatness

        def compute_contour_profile_flatness_radius_not_imposed(x):
            df = self.compute_2d_contour_scanned_profile_radius_not_imposed(
                angle=x[0],
                weight=x[1],
                shoot_at=x[2],
            )
            profile = df.query("x**2+y**2<=@self.fieldsize**2")["total_dose"]
            flatness = 1e2 * (_np.max(profile) - _np.min(profile)) / (_np.max(profile) + _np.min(profile))

            return flatness

        if self.shoot_on_aperture and not self.angle_imposed:
            optim = minimize(
                fun=compute_contour_profile_flatness,
                bounds=Bounds(
                    [0, 0],
                    [_np.pi, 10],
                ),
                x0=[_np.pi / 6, 0.5],
                method="trust-constr",
                options={"verbose": 1, "xtol": 1e-10, "maxiter": 1e5},
            )

        elif self.shoot_on_aperture and self.angle_imposed:
            optim = minimize(
                fun=compute_contour_profile_flatness_angle_imposed,
                bounds=Bounds(
                    [0],
                    [10],
                ),
                x0=0.5,
                method="trust-constr",
                options={"verbose": 1, "xtol": 1e-10, "maxiter": 1e5},
            )

        else:
            optim = minimize(
                fun=compute_contour_profile_flatness_radius_not_imposed,
                bounds=Bounds(
                    [0, 0, 0],
                    [_np.pi, 10, 2 * self.fieldsize],
                ),
                x0=[_np.pi / 6, 0.5, self.fieldsize],
                method="trust-constr",
                options={"verbose": 1, "xtol": 1e-10, "maxiter": 1e5},
            )

        return optim
