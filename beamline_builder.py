import pandas as pd

## TODO Move read_csv etc. from beamline to beamline_builder


class BeamlineBuilder:

    def replicate(a, n=2):
        # return n*a
        pass

    def join(a, b):
        # return a + b
        pass

    def __init__(self):
        self._beamline = []
        pass


    def add_bend(self, **kwargs):
        bend = {
            name: kwargs.get('name', 'bend'),
            angle: kwargs.get('angle', 0),
            K1: kwargs,get('K1', 0),
        }
        self._beamline.append(bend)
        return self

    def add_beamline(self):
        return self

    @property
    def line(self):
        return pd.DataFrame(self._beamline)




builder = BeamlineBUilder()
builder.add_bend(name="bend", a=1, b=2)
builder.line