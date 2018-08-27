import os
import pandas as pd
import numpy as np
import scipy.interpolate

STAR_COLUMNS = ['K', 'eS', 'nS', 'tS', 'csda', 'prange', 'factor']


class MaterialsDB:
    def __init__(self, p=os.path.dirname(__file__)):
        self.db_path = p
        pdg_data = self._pdg_read_data("data.csv")
        pstar_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(p, 'pstar')) if
                            os.path.isfile(os.path.join(p, 'pstar', f)) and os.path.splitext(f)[1] == '.txt']
        pstar = {
            m: self._star_read_data(m) for m in pstar_data_files
        }
        srim_data_files = [os.path.splitext(f)[0] for f in os.listdir(os.path.join(p, 'srim')) if
                           os.path.isfile(os.path.join(p, 'srim', f)) and os.path.splitext(f)[1] == '.txt']
        srim = {
            m: self._srim_read_data(m, pdg_data) for m in srim_data_files
        }
        projected_ranges_pstar = {
            k: scipy.interpolate.CubicSpline(np.log(v['K']), np.log(v['prange'])) for k, v in
            pstar.items()
        }
        projected_ranges_srim = {
            k: scipy.interpolate.CubicSpline(np.log(v['energy_scale']), np.log(v['projected_range_meter'])) for k, v in
            srim.items()
        }
        projected_ranges = {**projected_ranges_pstar, **projected_ranges_srim}
        csda_ranges_pstar = {
            k: scipy.interpolate.CubicSpline(np.log(v['K']), np.log(v['csda'])) for k, v in
            pstar.items()
        }
        csda_ranges = {**csda_ranges_pstar}
        self.__materials_db = {
            'projected_ranges': projected_ranges,
            'csda_ranges': csda_ranges,
            'pdg_data': pdg_data,
        }

    def _star_read_data(self, m):
        return pd.read_table(os.path.join(self.db_path, "pstar", f"{m}.txt"),
                             skiprows=7,
                             delimiter=' ',
                             names=STAR_COLUMNS,
                             index_col=False)

    def _srim_read_data(self, m, pdg_data):
        data = pd.read_csv(
            os.path.join(self.db_path, "srim", f"{m}.txt"),
            skiprows=23,
            skipfooter=20,
            engine='python',
            delim_whitespace=True,
            names=[
                'energy',
                'energy_unit',
                'stopping_elec',
                'stopping_nuclear',
                'projected_range',
                'projected_range_unit',
                'longitudinal_straggling',
                'longitudinal_straggling_unit',
                'lateral_straggling',
                'lateral_straggling_unit',
            ]
        )
        units = {
            'A': 1e-10,
            'um': 1e-6,
            'mm': 1e-3,
            'm': 1,
            'keV': 1e-3,
            'MeV': 1,
            'GeV': 1000,
        }

        data['projected_range'] = data['projected_range'].astype(float)
        data['projected_range_meter'] = data.apply(
            lambda e: units[e['projected_range_unit']] * pdg_data.at[m, 'rho'] * 100 * e['projected_range'], axis=1)
        data['energy'] = data['energy'].astype(float)
        data['energy_scale'] = data.apply(lambda e: units[e['energy_unit']] * e['energy'], axis=1)
        data['stopping_elec'] = data['stopping_elec'].astype(float)
        data['stopping_nuclear'] = data['stopping_nuclear'].astype(float)

        return data

    def _pdg_read_data(self, data):
        return pd.read_table(os.path.join(self.db_path, "pdg", data),
                             delimiter=',',
                             index_col='material'
                             )

    @property
    def db(self):
        return self.__materials_db

    @property
    def csda_ranges(self):
        return self.__materials_db['csda_ranges']

    @property
    def projected_ranges(self):
        return self.__materials_db['projected_ranges']

    def density(self, material):
        return self.__materials_db['pdg_data'].at[str(material), 'rho']

    def z(self, material):
        return self.__materials_db['pdg_data'].at[str(material), 'Z']

    def a(self, material):
        return self.__materials_db['pdg_data'].at[str(material), 'A']
