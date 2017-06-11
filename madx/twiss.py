import os
from madx.madx import read_madx_twiss
from madx.madx import read_ptc_twiss

def twiss(**kwargs):
    """Compute the Twiss parameters of the beamline."""
    line = kwargs.get('line', None)
    m = kwargs.get('madx', None)

    #self.__flag_ptc = kwargs.get('ptc', self.__flag_ptc)
    m.attach(line)
    m.beam()
    m.twiss(ptc=False, centre=True)
    errors = m.run(line.context).fatals
    line.__madx_input = m.input
    if len(errors) > 0:
        m.print_input()
        print(errors)
        raise "error"
        #raise BeamlineException("MAD-X ended with fatal error.")
    madx_twiss = read_madx_twiss(os.path.join(m.__path, 'twiss.outx'))
    line_with_twiss = madx_twiss.merge(line,
                                       left_index=True,
                                       right_index=True,
                                       how='outer',
                                       suffixes=('_TWISS', '')
                                       ).sort_values(by='S')
    return line_with_twiss
