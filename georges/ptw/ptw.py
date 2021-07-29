import numpy as _np
import pandas as pd
from numpy.polynomial import Polynomial
import scipy
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class BraggPeakAnalysis:

    def __init__(self, bp: pd.DataFrame, method: str, low_dose=1, high_dose=60):
        self.data = bp
        self.low_bound = low_dose
        self.high_bound = high_dose
        self.coefficients = None
        self.bp_interval = None
        self.method = method

    def compute_percentage(self, x):
        f = interp1d(self.data[self.data.columns[0]].values, self.data[self.data.columns[1]].values,
                     kind='linear', bounds_error=False)
        return f(x)

    def fitting_function(self, x):
        coeff = self.coefficients.convert().coef
        return coeff[0] + coeff[1] * x + coeff[2] * x ** 2 + coeff[3] * x ** 3 + coeff[4] * x ** 4 + \
               coeff[5] * x ** 5 + coeff[6] * x ** 6 + coeff[7] * x ** 7 + coeff[8] * x ** 8

    def fit_bp(self):
        """
        This function fits a Bragg peak using a 8th order polynom. Based on reference xxx.

        Args:
        Returns:

        """
        idx1 = _np.where(self.data[self.data.columns[1]] > self.high_bound)
        idx2 = _np.where(self.data[self.data.columns[1]] < self.low_bound)

        idx_min = idx1[0][0]
        idx_max = idx2[0][0]
        bp_interval = self.data[idx_min:idx_max]
        bp_interval.reset_index(drop=True, inplace=True)

        x = _np.float32(bp_interval[bp_interval.columns[0]])
        y = _np.float32(bp_interval[bp_interval.columns[1]])

        self.coefficients = Polynomial.fit(x=x,
                                           y=y,
                                           deg=8, full=False)

        self.bp_interval = bp_interval[bp_interval.columns[0]].to_numpy()

        return self.bp_interval, self.coefficients, self.fitting_function(self.bp_interval)

    def get_coefficients(self) -> _np.array:
        return self.coefficients.convert().coef

    def get_bp_interval(self) -> _np.array:
        return self.bp_interval

    def get_xrange(self, val):

        if self.method == 'polynom_fit':
            rval = (self.fit_bp()[1] - val).roots()
            rval = rval[rval.imag == 0].real
            var = rval[_np.where(_np.logical_and(rval >= self.bp_interval[0], rval <= self.bp_interval[
                -1]))]  # check if we are in the interval of the BP
            if var.size:
                return _np.max(var)
            return 0

        elif self.method == 'scipy.optimize':
            var = scipy.optimize.newton_krylov(lambda x: self.compute_percentage(x=x) - val,
                                               xin=self.data[self.data.columns[0]][
                                                   _np.abs(self.data[self.data.columns[1]] - val).idxmin()],
                                               maxiter=10000,
                                               f_tol=1e-12,
                                               x_tol=1e-13)
            return var
        else:
            raise Exception(
                "Please choose a valid computation method (must be either 'polynom_fit' or 'scipy.optimize')")

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
            _np.where(_np.logical_and(real_roots >= self.bp_interval[0], real_roots <= self.bp_interval[-1]))]
        fit_var = self.fitting_function(var)
        var = var[_np.where(fit_var == _np.max(fit_var))]  # take the value where the fit is maximum (global maximum)
        return var[0]

    def compute_distal_fall_off(self) -> float:
        return [self.get_r10() - self.get_r90(), self.get_r20() - self.get_r80()]

    def compute_dfo_90_10(self) -> float:
        return self.compute_distal_fall_off()[0]

    def compute_dfo_80_20(self) -> float:
        return self.compute_distal_fall_off()[1]


class SpreadOutBraggPeakAnalysis:

    def __init__(self, dose_data: pd.DataFrame, method: str, z_axis: _np.array, modul_type: str):
        self.dose_data = dose_data
        self.method = method
        self.z_axis = z_axis
        self.modul_type = 'Full'

    def max_indexes(self) -> _np.array:
        max_indexes = _np.array(list(map(int, self.dose_data.idxmax(axis=1).values.tolist())))

        return max_indexes

    def sobp_data(self) -> _np.array:
        sobp_data = _np.zeros((self.dose_data.shape[0], self.dose_data.shape[0]))

        for i in range(self.dose_data.shape[0]):
            sobp_data[i] = self.dose_data.iloc[:, self.max_indexes()[i]]

        return sobp_data

    def compute_weights(self) -> _np.array:
        goal_dose_values = _np.full(self.dose_data.shape[0], 1)  # Here we want 1 Gy

        if self.method == 'np.linalg.solve':

            return _np.linalg.solve(self.sobp_data(), goal_dose_values)

        elif self.method == 'scipy.optimize':

            a = self.sobp_data()
            b = goal_dose_values
            n = len(b)
            fun = lambda x: _np.linalg.norm(_np.dot(a, x) - b)
            sol = minimize(fun,
                           _np.zeros(n),
                           method='L-BFGS-B',
                           bounds=[(0., None) for x in range(n)])

            return sol['x']

    def compute_sobp_profile(self) -> _np.array:

        weighted_dose_data = _np.matmul(self.compute_weights().reshape(1, self.dose_data.shape[0]),
                                        self.dose_data.values)

        return weighted_dose_data.reshape(self.dose_data.shape[1])

    def view_sobp(self, with_pristine_peaks=False):

        def plot(x, y):
            plt.plot(x,
                     y,
                     linestyle='dashed',
                     linewidth=1.2,
                     color='r',
                     marker='*',
                     markersize=2,
                     )

        if with_pristine_peaks:
            weighted_dose = self.dose_data.multiply(self.compute_weights(), axis=0)
            weighted_dose.apply(lambda e: plot(x=self.z_axis,
                                               y=1e2 * e.values),
                                axis='columns')

        plt.plot(self.z_axis,
                 100 * self.compute_sobp_profile(),
                 linestyle='dashed',
                 linewidth=2,
                 color='k',
                 marker='*',
                 markersize=3,
                 label='SOBP')

        plt.xlabel('Depth (mm)')
        plt.ylabel('Normalized dose (\\%)')
        plt.legend(loc='center left', ncol=1)
        plt.xlim(0, _np.max(self.z_axis))
        plt.ylim(0, _np.max(100 * self.compute_sobp_profile()))
        plt.grid()

    def recompute_sobp_for_another_range(self, bp_to_leave):

        if bp_to_leave == 0:
            new_weighted_dose = _np.matmul(
                self.compute_weights().reshape(1, self.dose_data.shape[0]),
                self.dose_data.values[0:, :])
        else:
            new_weighted_dose = _np.matmul(self.compute_weights()[:-bp_to_leave].reshape(1,
                                                                                         self.dose_data.shape[
                                                                                             0] - bp_to_leave),
                                           self.dose_data.values[bp_to_leave:, :])

        sobp_array = 1e2 * new_weighted_dose.reshape(self.dose_data.shape[1])

        return sobp_array

    def recompute_sobp_for_another_width(self, bp_to_leave):

        weights = _np.flip(_np.flip(self.compute_weights())[bp_to_leave:])
        weighted_dose_data = _np.matmul(weights.reshape(1, self.dose_data.shape[0] - bp_to_leave),
                                        self.dose_data.values[bp_to_leave:, :])

        sobp_data = 1e2 * weighted_dose_data.reshape(self.dose_data.shape[1])

        return sobp_data

    def compute_ranges_and_flatness(self):

        sobp_array = 1e2 * self.compute_sobp_profile()
        idxmax = _np.where(sobp_array < 1)[0]
        idx = 0

        for i in range(_np.flip(sobp_array[0:idxmax[0]]).shape[0]):
            if _np.flip(sobp_array[0:idxmax[0]])[i + 1] > _np.flip(sobp_array[0:idxmax[0]])[i]:
                idx += 1
            else:
                break

        flipped_data = _np.flip(sobp_array[idxmax[0] - idx:idxmax[0]])

        while _np.max(flipped_data) < 98:
            idx += 1
            flipped_data = _np.flip(sobp_array[idxmax[0] - idx:idxmax[0]])

        def compute_range(data, r, z_axis):

            df_idx = (pd.DataFrame(_np.abs(data - r))).idxmin()
            xin = _np.flip(z_axis[idxmax[0] - idx:idxmax[0]])[df_idx][0]

            function = interp1d(_np.flip(z_axis[idxmax[0] - idx:idxmax[0]]),
                                _np.flip(sobp_array[idxmax[0] - idx:idxmax[0]]).reshape(idx),
                                kind='linear',
                                bounds_error=False)
            res = scipy.optimize.newton_krylov(lambda x: function(x=x) - r,
                                                xin=xin,
                                                maxiter=1000000,
                                                f_tol=1e-1,
                                                x_tol=1e-1)
            return res

        r_10 = compute_range(flipped_data, 10, self.z_axis)
        r_90 = compute_range(flipped_data, 90, self.z_axis)
        r_98 = compute_range(flipped_data, 98, self.z_axis)

        if self.modul_type == 'Full':
            max_dose = _np.max(sobp_array[_np.where(self.z_axis < r_98)[0]])
            min_dose = _np.min(sobp_array[_np.where(self.z_axis < r_98)[0]])
            flatness = 100 * ((max_dose - min_dose) / (max_dose + min_dose))

            return r_10, r_90, r_98, flatness

    def get_sobp_r10(self):
        return self.compute_ranges_and_flatness()[0]

    def get_sobp_r90(self):
        return self.compute_ranges_and_flatness()[1]

    def get_sobp_r98(self):
        return self.compute_ranges_and_flatness()[2]

    def get_sobp_flatness(self):
        return self.compute_ranges_and_flatness()[3]


class LateralProfileAnalysis:

    def __init__(self, dose_profile: _np.array, positions: _np.array):
        self.dose_profile = dose_profile
        self.positions = positions

    def set_data(self):
        idxmax = _np.where(self.positions > 0)[0][0]

        positions_left = self.positions[0:idxmax]
        positions_right = self.positions[idxmax:]
        dose_left = self.dose_profile[0:idxmax]
        dose_right = self.dose_profile[idxmax:]

        return idxmax, positions_left, dose_left, positions_right, dose_right

    def define_f(self):
        f_left = interp1d(self.set_data()[1],
                          self.set_data()[2],
                          kind='cubic', bounds_error=False)
        f_right = interp1d(self.set_data()[3],
                           self.set_data()[4],
                           kind='cubic', bounds_error=False)
        return f_left, f_right

    def get_position_left(self, percentage):

        position = scipy.optimize.newton_krylov(lambda x: self.define_f()[0](x=x) - percentage,
                                                xin=self.set_data()[1][
                                                     _np.abs(pd.DataFrame(self.set_data()[2] - percentage)).idxmin()],
                                                maxiter=100000,
                                                f_tol=1e-2,
                                                x_tol=1e-2)
        return position[0]

    def get_position_right(self, percentage):

        position = scipy.optimize.newton_krylov(lambda x: self.define_f()[1](x=x) - percentage,
                                                xin=self.set_data()[3][
                                                      _np.abs(pd.DataFrame(self.set_data()[4] - percentage)).idxmin()],
                                                maxiter=100000,
                                                f_tol=1e-2,
                                                x_tol=1e-2)
        return position[0]

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
        return 1e2 * (self.get_ur_max_dose() - self.get_ur_min_dose()) / (
                self.get_ur_max_dose() + self.get_ur_min_dose())

    def get_penumbra_left(self):
        return self.get_p_80_left() - self.get_p_20_left()

    def get_penumbra_right(self):
        return self.get_p_20_right() - self.get_p_80_right()

    def get_penumbra(self):
        return _np.mean([self.get_penumbra_left(), self.get_penumbra_right()])