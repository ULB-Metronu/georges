import zgoubidoo
from .converters import Converters
from .. import Beamline


def convert(beamline: Beamline):
    flag = False
    tmp = list()
    for e in Beamline(
            beamline.line.query("TYPE != 'MARKER'")
        ).add_drifts(using_collimators=True).line.apply(
            lambda _: getattr(Converters, e.TYPE.lower(), None)(_),
            axis=1
        ).dropna().values.tolist():
        if isinstance(e, zgoubidoo.commands.Dipole):
            if e.AT < 0:
                if flag is False:
                    b.append(Ymy())
                    flag = True
                e.AT *= -1
        tmp.append(e)

    return zgoubidoo.Input(beamline.name, tmp)
