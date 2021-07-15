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

        # see https://stackoverflow.com/questions/13195832/numpy-ma-polyfit-function-for-masked-arrays-crashes-on-integer-input
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
            var = rval[_np.where(_np.logical_and(rval >= self.bp_interval[0], rval <= self.bp_interval[-1]))] # check if we are in the interval of the BP
            if var.size:
                return _np.max(var)
            return 0

        elif self.method == 'scipy.optimize':
            var = scipy.optimize.newton_krylov(lambda x: self.compute_percentage(x=x) - val, \
                                               xin=self.data[self.data.columns[0]][_np.abs(self.data[self.data.columns[1]]-val).idxmin()],
                                               maxiter=10000,
                                               f_tol=1e-12,
                                               x_tol=1e-13)
            return var
        else:
            raise Exception("Please choose a valid computation method (must be either 'polynom_fit' or 'scipy.optimize')")


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
        var = real_roots[_np.where(_np.logical_and(real_roots >= self.bp_interval[0], real_roots <= self.bp_interval[-1]))]
        fit_var = self.fitting_function(var)
        var = var[_np.where(fit_var == _np.max(fit_var))] # take the value where the fit is maximum (global maximum)
        return var[0]

    def compute_distal_fall_off(self) -> float:
        return [self.get_r10() - self.get_r90(), self.get_r20() - self.get_r80()]

    def compute_dfo_90_10(self) -> float:
        return self.compute_distal_fall_off()[0]

    def compute_dfo_80_20(self) -> float:
        return self.compute_distal_fall_off()[1]


class SpreadOutBraggPeakAnalysis:

    def __init__(self, dose_data: pd.DataFrame, method: str, z_axis: _np.array):
        self.dose_data = dose_data
        self.method = method
        self.z_axis = z_axis

    def max_indexes(self) -> _np.array:
        max_indexes = _np.array(list(map(int, self.dose_data.idxmax(axis=1).values.tolist())))

        return max_indexes

    def sobp_data(self) -> _np.array:
        sobp_data = _np.zeros((self.dose_data.shape[0],self.dose_data.shape[0]))

        for i in range(self.dose_data.shape[0]):
            sobp_data[i] = self.dose_data.iloc[:,self.max_indexes()[i]]

        return sobp_data

    def compute_weights(self) -> _np.array:
        goal_dose_values = _np.full(self.dose_data.shape[0], 1) # Here we want 1 Gy

        if self.method == 'np.linalg.solve':
            return _np.linalg.solve(self.sobp_data(), goal_dose_values)

        elif self.method == 'scipy.optimize':

            A = self.sobp_data()
            b = goal_dose_values
            n = len(b)
            fun = lambda x: _np.linalg.norm(_np.dot(A,x)-b)
            sol = minimize(fun,
                           _np.zeros(n),
                           method='L-BFGS-B',
                           bounds=[(0.,None) for x in range(n)])

            return sol['x']

    def compute_sobp_profile(self) -> _np.array:

        weighted_dose_data = _np.matmul(self.compute_weights().reshape(1,self.dose_data.shape[0]),
                                       self.dose_data.values)

        return weighted_dose_data.reshape(self.dose_data.shape[1])

    def view_sobp(self, with_pristine_peaks = False):
        if with_pristine_peaks == True:
            for i in range(self.dose_data.shape[0]-1):
                plt.plot(self.z_axis,
                         100*self.compute_weights()[i]*self.dose_data.values[i,:],
                         linestyle='dashed',
                         linewidth=1.2,
                         color='r',
                         marker='*',
                         markersize=2,
                         )
            plt.plot(self.z_axis,
                     100*self.compute_weights()[self.dose_data.shape[0]-1]*self.dose_data.values[self.dose_data.shape[0]-1,:],
                     linestyle='dashed',
                     linewidth=1.2,
                     color='r',
                     marker='*',
                     markersize=2,
                     label='Pristine peaks'
                     )

        plt.plot(self.z_axis,
                 100*self.compute_sobp_profile(),
                 linestyle='dashed',
                 linewidth=2,
                 color='k',
                 marker='*',
                 markersize=3,
                 label='SOBP')
        plt.xlabel('Depth (mm)')
        plt.ylabel('Normalized dose (\\%)')
        plt.legend(loc='center left',ncol=1)
        plt.xlim(0,_np.max(self.z_axis))
        plt.ylim(0,_np.max(100*self.compute_sobp_profile()))
        plt.grid()

    
class LateralProfileAnalysis:
    
    def __init__(self, dose_profile: _np.array, positions: _np.array):
        self.dose_profile = dose_profile
        self.positions = positions
        
    def set_data(self):
        
        idxmax = _np.where(self.positions>0)[0][0]

        positions_left = self.positions[0:idxmax]
        dose_left = self.dose_profile[0:idxmax]

        positions_right = self.positions[idxmax:]
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

    
    def get_p_20_left(self):
        p_20_left = scipy.optimize.newton_krylov(lambda x: self.define_f()[0](x=x) - 20,
                                xin=self.set_data()[1][_np.abs(pd.DataFrame(self.set_data()[2]-20)).idxmin()],
                                maxiter=100000,
                                f_tol=1e-2,
                                x_tol=1e-2)
        return p_20_left[0]
   
    def get_p_50_left(self):
        p_50_left = scipy.optimize.newton_krylov(lambda x: self.define_f()[0](x=x) - 50,
                                        xin=self.set_data()[1][_np.abs(pd.DataFrame(self.set_data()[2]-50)).idxmin()],
                                        maxiter=100000,
                                        f_tol=1e-2,
                                        x_tol=1e-2)
        return p_50_left[0]
    
    def get_p_80_left(self):
        p_80_left = scipy.optimize.newton_krylov(lambda x: self.define_f()[0](x=x) - 80,
                                    xin=self.set_data()[1][_np.abs(pd.DataFrame(self.set_data()[2]-80)).idxmin()],
                                    maxiter=100000,
                                    f_tol=1e-2,
                                    x_tol=1e-2)
        return p_80_left[0]
    
    def get_p_20_right(self):
        p_20_right = scipy.optimize.newton_krylov(lambda x: self.define_f()[1](x=x) - 20,
                                    xin=self.set_data()[3][_np.abs(pd.DataFrame(self.set_data()[4]-20)).idxmin()],
                                    maxiter=100000,
                                    f_tol=1e-2,
                                    x_tol=1e-2)
        return p_20_right[0]
    
    def get_p_50_right(self):
        p_50_right = scipy.optimize.newton_krylov(lambda x: self.define_f()[1](x=x) - 50,
                                    xin=self.set_data()[3][_np.abs(pd.DataFrame(self.set_data()[4]-50)).idxmin()],
                                    maxiter=100000,
                                    f_tol=1e-2,
                                    x_tol=1e-2)
        return p_50_right[0]
    
    def get_p_80_right(self):
        p_80_right = scipy.optimize.newton_krylov(lambda x: self.define_f()[1](x=x) - 80,
                                    xin=self.set_data()[3][_np.abs(pd.DataFrame(self.set_data()[4]-80)).idxmin()],
                                    maxiter=100000,
                                    f_tol=1e-2,
                                    x_tol=1e-2)
        return p_80_right[0]
    
    def get_field_size(self):
        return self.get_p_50_right() - self.get_p_50_left()
    
    def get_ur_left(self):
        return self.get_p_50_left() + 3
    
    def get_ur_right(self):
        return self.get_p_50_right() - 3    
    
    def get_ur_size(self):
        return self.ur_right() - self.get_ur_left()


    def get_ur_min_dose(self):
        ur_idx_min = _np.where(self.get_ur_left()<self.positions)[0][0]
        ur_idx_max = _np.where(self.positions>self.get_ur_right())[0][0]
        
        return _np.min(self.dose_profile[ur_idx_min:ur_idx_max])

    def get_ur_max_dose(self):
        ur_idx_min = _np.where(self.get_ur_left()<self.positions)[0][0]
        ur_idx_max = _np.where(self.positions>self.get_ur_right())[0][0]
        
        return _np.max(self.dose_profile[ur_idx_min:ur_idx_max])

    def get_ur_flatness(self):
        return 1e2 * (self.get_ur_max_dose() - self.get_ur_min_dose()) / (self.get_ur_max_dose() + self.get_ur_min_dose())

    def get_penumbra_left(self):
        return self.get_p_80_left() - self.get_p_20_left()

    def get_penumbra_right(self):
        return self.get_p_20_right() - self.get_p_80_right()

    def get_penumbra(self):
        return _np.mean([self.get_penumbra_left(), self.get_penumbra_right()])