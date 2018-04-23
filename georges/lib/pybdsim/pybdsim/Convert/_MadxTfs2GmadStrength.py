from .. import Builder as _Builder
from .. import _General as _General
from . import _MadxTfs2Gmad
import pymadx as _pymadx

_ElementModifier = _Builder.ElementModifier
_lFake           = _MadxTfs2Gmad._lFake

def MadxTfs2GmadStrength(input, outputfilename, existingmachine=None, verbose=False, flipmagnets=False):
    """
    Use a MADX Tfs file containing full twiss information to generate strength (only) file
    to be used with an existing lattice.

    """
    # ensure it's tfs instance and if not open the filepath provided
    madx = _pymadx.CheckItsTfs(input)

    # Zero any missing required columns
    ZeroMissingRequiredColumns(madx)

    if existingmachine == None:
        existingnames = madx.GetColumn('NAME')
    else:
        existingnames = existingmachine.keys() #dictionary

    if verbose:
        madx.ReportPopulations()

    #existing machine can be a pybdsim.Builder.Machine instance or a list or a dictionary instance
    #but should contain the desired keys that strengths should be generated for from the Tfs file.

    newStrengths = []
    
    for item in madx:
        name  = item['NAME']
        rname = _General.PrepareReducedName(name)
        
        # generate elementmodifier with approprate name to match one
        # already used in existing machine
        if name in existingnames:
            a = _GenerateElementModifier(item, name, verbose, flipmagnets)
        elif rname in existingnames:
            a = _GenerateElementModifier(item, rname, verbose, flipmagnets)
        else:
            a = None
        if verbose:
            print(a)
        if a != None:
            newStrengths.append(a)

    #write output
    if not outputfilename.endswith('.gmad'):
        outputfilename += '.gmad'
    f = open(outputfilename, 'w')
    for strength in newStrengths:
        f.write(str(strength))
    f.close()


def _GenerateElementModifier(madxitem, nameToUse, verbose=False, flipmagnets=False):
    """
    Generate an Builder.ElementModifier instance based on the particular
    element / magnet type.  

    Takes second argument of nameToUse to match whatever the name as been
    changed to.

    This function doesn't do any key checking for the dictionary as that should
    be done by MadxTfs2GmadStrength
    """
    item = madxitem
    category = item['KEYWORD']

    l      = item['L']
    factor = -1 if flipmagnets else 1  #flipping magnets
    
    a = None
    if category == 'SOLENOID':
        newk0 = item['K0L'] / l  * factor
        a = _ElementModifier(nameToUse, ks=newk0)
    elif category == 'QUADRUPOLE':
        newk1 = item['K1L'] / l * factor
        a = _ElementModifier(nameToUse, k1=newk1)
    elif category == 'SEXTUPOLE':
        newk2 = item['K2L'] / l * factor
        a = _ElementModifier(nameToUse, k2=newk2)
    elif category == 'OCTUPOLE':
        newk3 = item['K3L'] / l * factor
        a = _ElementModifier(nameToUse, k3=newk3)
    elif category == 'DECAPOLE':
        newk4 = item['K4L'] / l * factor
        a = _ElementModifier(nameToUse, k4=newk4)
    elif category == 'HKICKER':
        a = _ElementModifier(nameToUse, angle=item['HKICK']*factor)
    elif category == 'VKICKER':
        a = _ElementModifier(nameToUse, angle=item['VKICK']*factor)
    elif category == 'TKICKER':
        if verbose:
            print('WARNING - TKICKER not implemented yet')
    elif category == 'RFCAVITY':
        if verbose:
            print('WARNING - RFCAVITY not implemented yet')
    elif category == 'MULTIPOLE':
        k1  = item['K1L']  / _lFake * factor
        k2  = item['K2L']  / _lFake * factor
        k3  = item['K3L']  / _lFake * factor
        k4  = item['K4L']  / _lFake * factor
        k5  = item['K5L']  / _lFake * factor
        k6  = item['K6L']  / _lFake * factor
        k1s = item['K1SL'] / _lFake * factor
        k2s = item['K2SL'] / _lFake * factor
        k3s = item['K3SL'] / _lFake * factor
        k4s = item['K4SL'] / _lFake * factor
        k5s = item['K5SL'] / _lFake * factor
        k6s = item['K6SL'] / _lFake * factor
        a = _ElementModifier(nameToUse, True, knl=(k1,k2,k3,k4,k5,k6), ksl=(k1s,k2s,k3s,k4s,k5s,k6s))
    else:
        pass # just keep a = None and return that

    if verbose:
        print('New Strength: ', a)
        if a == None:
            print('Unsupported type: ',item['KEYWORD'])
    
    return a
