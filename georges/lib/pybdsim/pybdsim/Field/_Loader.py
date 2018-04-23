import numpy as _np
import tarfile as _tarfile

def Load(filename, debug=False):
    """
    Load a BDSIM field format file into a numpy array. Can either
    be a regular ascii text file or can be a compressed file ending
    in ".tar.gz".

    returns a numpy array with the corresponding number of dimensions
    and the dimension has the coordaintes and fx,fy,fz.
    """
    if (filename.endswith('.tar.gz')):
        print('Field Loader> loading compressed file ' + filename)
        tar = _tarfile.open(filename,'r')
        f = tar.extractfile(tar.firstmember)
    else:
        print('Field Loader> loading file ' + filename)
        f = open(filename)

    intoData = False
    header   = {}
    columns  = []
    data     = []
    
    for line in f:
        if intoData:
            data.append(line.strip().split())

        elif '>' in line:
            d = line.strip().split('>')
            k = d[0].strip()
            v = float(d[1].strip())
            header[k] = v

        elif '!' in line:
            columns = line.strip('!').strip().split()
            intoData = True

    f.close()

    data = _np.array(data, dtype=float)

    nDim = len(columns) - 3
    if (nDim < 1 or nDim > 4):
        if debug:
            print('Invalid number of columns')
            print(columns)
        return

    required = ['nx','ny','nz','nt']

    requiredKeys    = required[:nDim]
    requiredKeysSet = set(requiredKeys)
    if not requiredKeysSet.issubset(header.keys()):
        print('missing keys from header!')
        if debug:
            print(header)
        return
    else:
        dims = [int(header[k]) for k in requiredKeys[::-1]]
        dims.append(len(columns))
        if debug:
            print(dims)
            print(nDim)
            print(_np.shape(data))
        data = data.reshape(*dims)

    return data
