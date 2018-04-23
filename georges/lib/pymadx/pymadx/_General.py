# pymadx._General - general python scripts / tools
# Version 1.0
# L. Nevay, S.T.Boogert
# laurie.nevay@rhul.ac.uk

"""
General utilities for day to day housekeeping
"""

import os

def CheckFileExists(filename):
    i = 1
    parts = filename.split('.')
    basefilename = parts[0]
    if len(parts) > 1:
        extension = '.' + parts[1]
    else:
        extension = ''
    while os.path.exists(filename) :
        filename = basefilename+str(i)+extension
        i = i + 1
    return filename

def Chunks(l, n):
    """ Yield successive n-sized chunks from l.    """
    return [l[i:i+n] for i in range(0,len(l),n)]

def NearestEvenInteger(number):
    number = int(number)
    return number + number%2

def Cast(string):
    """
    Cast(string)
    
    tries to cast to a (python)float and if it doesn't work, 
    returns a string

    """
    try:
        return float(string)
    except ValueError:
        return string

def IsFloat(stringtotest):
    try:
        float(stringtotest)
        return True
    except ValueError:
        return False

def IndexOfElement(tfsinstance,markername):
    t = tfsinstance
    names = list(t.data['NAME'])
    try:
        i = names.index(markername)
    except ValueError:
        i = 0
        print('Unknown element name')
    return i

def GetSixTrackAperType(aper1,aper2,aper3,aper4):
    if aper1 == 0 and aper2 == 0 and aper3 == 0 and aper4== 0:
        return ''
    elif aper1 == aper3 and aper2 == aper4:
        return 'ELLIPSE'
    elif aper1 == aper3 and aper2 < aper4:
        return 'LHCSCREEN'
    elif aper1 < aper3 and aper2 == aper4:
        return 'LHCSCREEN'
    elif aper1 == 0 and aper2 == 0:
        return 'RACETRACK'
    elif aper3 == 0:
        return 'RECTANGLE'
    else:
        s = "WARNING: The given aperture is not classified among the known types\n"
        s += "A1 = " + str(aper1) + ", A2 = " +  str(aper2) + ", A3 = "
        s += str(aper3) + ", A4 = " + str(aper4)
        raise AttributeError(s)
