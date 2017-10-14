import pandas as pd
from .physics import *


def read_range_file(self, file):
    data = pd.read_csv(file)
    self.__strength = data


class Context:

    def __init__(self):
        self.__context = {
            'PARTICLE': 'PROTON',

        }
        self.__strength = None



    def __getitem__(self, key):
        return self.__context[key]

