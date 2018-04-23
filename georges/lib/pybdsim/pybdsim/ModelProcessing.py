# pybdsim.ModelProcessing - tools to process existing BDSIM models
# Version 1.0
# L. Nevay
# laurie.nevay@rhul.ac.uk

"""
ModelProcessing

Tools to process existing BDSIM models and generate other
versions of them.

"""

from . import _General
from . import Gmad as _Gmad
from . import Builder as _Builder

import time as _time

def GenerateFullListOfSamplers(inputfile, outputfile):
    """
    inputfile - path to main gmad input file

    This will parse the input using the compiled BDSIM
    parser (GMAD), iterate over all the beamline elements
    and generate a sampler for every elements.  Ignores
    samplers, but may include already defined ones in your
    own input.

    """
    lattice  = _Gmad.Lattice(inputfile)
    samplers = []
    typestoignore = [
        'None',
        'Sampler',
        'CSampler',
        'Line',
        'Reversed Line',
        'Material',
        'Atom',
        'Sequence',
        'Tunnel',
        'Teleporter',
        'Terminator',
        'Transform3D'
        ]
    for e in lattice:
        if e['Type'] not in typestoignore:
           samplers.append(_Builder.Sampler(e['Name']))

    _WriteSamplerToGmadFile(samplers,outputfile)

def _WriteSamplerToGmadFile(samplerlist, outputfile):
    f = open(outputfile, "w")
    for s in samplerlist:
        f.write(s.__repr__())
    f.close()
                           

def WrapLatticeAboutItem(maingmadfile, itemname, outputfilename):
    elementsperline = 100
    
    a = _Gmad.Lattice(maingmadfile)

    seq = a.sequence

    # remove 'lattice' from list
    try:
        seq.remove('lattice')
    except ValueError:
        pass

    # remove samplers from list
    seq = [name for name in seq if 'ampler' not in name]

    try:
        ind = seq.index(itemname)
    except ValueError:
        print('Error - ',itemname,'is not in the supplied lattice')
        return

    newseq = seq[ind:] + seq[:ind]

    # write lattice definition
    f = open(outputfilename,'w')
    timestring = '! ' + _time.strftime("%a, %d %b %Y %H:%M:%S +0000", _time.gmtime()) + '\n'
    f.write(timestring)
    f.write('! pybdsim.ModelProcessing Wrapped Lattice about '+str(itemname) +'\n')
    f.write('! LATTICE SEQUENCE DEFINITION\n\n')
    sequencechunks = _General.Chunks(newseq,elementsperline)
    linelist = []
    ti = 0
    for line in sequencechunks:
        f.write('l'+str(ti)+': line = ('+', '.join(line)+');\n')
        linelist.append('l'+str(ti))
        ti += 1
    # need to define the period before making sampler planes
    f.write('lattice: line = ('+', '.join(linelist)+');\n')
    f.write('use, lattice;\n')
    f.close()
