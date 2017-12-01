import os
import pandas as pd
import numpy as np
import scipy.interpolate

STAR_COLUMNS = ['K', 'eS', 'nS', 'tS', 'csda', 'prange', 'factor']


class MaterialsDB:
    def __init__(self, p):
        self.db_path = p
        data_files = [os.path.splitext(f)[0] for f in os.listdir(p) if
                      os.path.isfile(os.path.join(p, f)) and os.path.splitext(f)[1] == '.txt']
        pstar = {m: self._star_read_data(m) for m in data_files}
        projected_ranges = {k: scipy.interpolate.CubicSpline(np.log(v['K']), np.log(v['prange'])) for k, v in
                            pstar.items()}
        csda_ranges = {k: scipy.interpolate.CubicSpline(np.log(v['K']), np.log(v['csda'])) for k, v in pstar.items()}
        self.__materials_db = {
            'projected_ranges': projected_ranges,
            'csda_ranges': csda_ranges,
        }

    def _star_read_data(self, m):
        data = pd.read_table(os.path.join(self.db_path, f"{m}.txt"),
                             skiprows=7,
                             delimiter=' ',
                             names=STAR_COLUMNS,
                             index_col=False)

        # data.index.name = 'K'
        return data

    @property
    def csda_ranges(self):
        return self.__materials_db['csda_ranges']

    @property
    def projected_ranges(self):
        return self.__materials_db['projected_ranges']

    @staticmethod
    def density(material):
        density = {
            'water': 1,
            'graphite': 1.7,
        }
        return density[material]
