import numpy as _np


class Field(object):
    """
    Base class used for common writing procedures for BDSIM field format.

    This does not support arbitrary loop ordering - only the originally intended
    xyzt.
    """
    def __init__(self, array=_np.array([]), columns=[], flip=True, doublePrecision=False):
        self.data            = array
        self.columns         = columns
        self.header          = {}
        self.flip            = flip
        self.doublePrecision = doublePrecision       

    def Write(self, fileName):
        f = open(fileName, 'w')
        for key,value in self.header.iteritems():
            f.write(str(key)+'> '+ str(value) + '\n')

        if self.doublePrecision:
            colStrings = map(lambda s: '%23s' % s, self.columns)
        else:
            colStrings = map(lambda s: '%14s' % s, self.columns)
        colStrings[0] = colStrings[0].strip() # don't pad the first column title
        # a '!' denotes the column header line
        f.write('! '+ '\t'.join(colStrings)+'\n')
        
        # flatten all but last dimension - 3 field components
        nvalues = _np.shape(self.data)[-1] # number of values in last dimension

        if self.flip:
            # [x,y,z,t,values] -> [t,z,y,x,values] for 4D
            # [x,y,z,values]   -> [z,y,x,values]   for 3D
            # [x,y,values]     -> [y,x,values]     for 2D
            # [x,values]       -> [x,values]       for 1D
            if (self.data.ndim == 2):
                pass # do nothin for 1D
            inds = range(self.data.ndim)       # indices for dimension [0,1,2] etc
            # keep the last value the same but reverse all indices before then
            inds[:(self.data.ndim - 1)] = reversed(inds[:(self.data.ndim - 1)])
            datal = _np.transpose(self.data, inds)
        else:
            datal = self.data

        datal = datal.reshape(-1,nvalues)
        for value in datal:
            if self.doublePrecision:
                strings   = map(lambda x: '%.16E' % x, value)
                stringsFW = map(lambda s: '%23s' % s,  strings)
            else:
                strings   = map(lambda x: '%.8E' % x, value)
                stringsFW = map(lambda s: '%14s' % s, strings)
            f.write('\t'.join(stringsFW) + '\n')

        f.close()


class Field1D(Field):
    """
    Utility class to write a 1D field map array to BDSIM field format.

    The array supplied should be 2 dimensional. Dimensions are:
    (x,value) where value has 4 elements [x,fx,fy,fz]. So a 120 long
    array would have np.shape of (120,4).
    
    This can be used for both electric and magnetic fields.

    Example::
    
    >>> a = Field1D(data)
    >>> a.Write('outputFileName.dat')

    """
    def __init__(self, data, doublePrecision=False):
        columns = ['X','Fx','Fy','Fz']
        super(Field1D, self).__init__(data,columns,doublePrecision=doublePrecision)
        self.header['xmin'] = _np.min(self.data[:,0])
        self.header['xmax'] = _np.max(self.data[:,0])
        self.header['nx']   = _np.shape(self.data)[0]

class Field2D(Field):
    """
    Utility class to write a 2D field map array to BDSIM field format.

    The array supplied should be 3 dimensional. Dimensions are:
    (x,y,value) where value has 5 elements [x,y,fx,fy,fz].  So a 100x50 (x,y)
    grid would have np.shape of (100,50,5).

    Example::
    
    >>> a = Field2D(data) # data is a prepared array
    >>> a.Write('outputFileName.dat')

    The 'flip' boolean allows an array with (y,x,value) dimension order
    to be written as (x,y,value).

    The 'doublePrecision' boolean controls whether the field and spatial
    values are written to 16 s.f. (True) or 8 s.f. (False - default).

    """
    def __init__(self, data, flip=True, doublePrecision=False):
        columns = ['X','Y','Fx','Fy','Fz']
        super(Field2D, self).__init__(data,columns,flip,doublePrecision)
        self.header['xmin'] = _np.min(self.data[:,:,0])
        self.header['xmax'] = _np.max(self.data[:,:,0])
        self.header['nx']   = _np.shape(self.data)[1]
        self.header['ymin'] = _np.min(self.data[:,:,1])
        self.header['ymax'] = _np.max(self.data[:,:,1])
        self.header['ny']   = _np.shape(self.data)[0]

class Field3D(Field):
    """
    Utility class to write a 3D field map array to BDSIM field format.

    The array supplied should be 4 dimensional. Dimensions are:
    (x,y,z,value) where value has 6 elements [x,y,z,fx,fy,fz].  So a 100x50x30 
    (x,y,z) grid would have np.shape of (100,50,30,6).
    
    Example::
    
    >>> a = Field3D(data) # data is a prepared array
    >>> a.Write('outputFileName.dat')

    The 'flip' boolean allows an array with (z,y,x,value) dimension order to
    be written as (x,y,z,value).

    The 'doublePrecision' boolean controls whether the field and spatial
    values are written to 16 s.f. (True) or 8 s.f. (False - default).

    """
    def __init__(self, data, flip=True, doublePrecision=False):
        columns = ['X','Y','Z','Fx','Fy','Fz']
        super(Field3D, self).__init__(data,columns,flip,doublePrecision)
        self.header['xmin'] = _np.min(self.data[:,:,:,0])
        self.header['xmax'] = _np.max(self.data[:,:,:,0])
        self.header['nx']   = _np.shape(self.data)[2]
        self.header['ymin'] = _np.min(self.data[:,:,:,1])
        self.header['ymax'] = _np.max(self.data[:,:,:,1])
        self.header['ny']   = _np.shape(self.data)[1]
        self.header['zmin'] = _np.min(self.data[:,:,:,2])
        self.header['zmax'] = _np.max(self.data[:,:,:,2])
        self.header['nz']   = _np.shape(self.data)[0]

class Field4D(Field):
    """
    Utility class to write a 4D field map array to BDSIM field format.

    The array supplied should be 5 dimensional. Dimensions are:
    (t,y,z,x,value) where value has 7 elements [x,y,z,t,fx,fy,fz]. So a 100x50x30x10
    (x,y,z,t) grid would have np.shape of (10,30,50,100,7).
    
    Example::
    
    >>> a = Field4D(data) # data is a prepared array
    >>> a.Write('outputFileName.dat')

    The 'flip' boolean allows an array with (t,z,y,x,value) dimension order to
    be written as (x,y,z,t,value).

    The 'doublePrecision' boolean controls whether the field and spatial
    values are written to 16 s.f. (True) or 8 s.f. (False - default).

    """
    def __init__(self, data, flip=True, doublePrecision=False):
        columns = ['X','Y','Z','T','Fx','Fy','Fz']
        super(Field4D, self).__init__(data,columns,flip,doublePrecision)
        self.header['xmin'] = _np.min(self.data[:,:,:,:,0])
        self.header['xmax'] = _np.max(self.data[:,:,:,:,0])
        self.header['nx']   = _np.shape(self.data)[3]
        self.header['ymin'] = _np.min(self.data[:,:,:,:,1])
        self.header['ymax'] = _np.max(self.data[:,:,:,:,1])
        self.header['ny']   = _np.shape(self.data)[2]
        self.header['zmin'] = _np.min(self.data[:,:,:,:,2])
        self.header['zmax'] = _np.max(self.data[:,:,:,:,2])
        self.header['nz']   = _np.shape(self.data)[1]
        self.header['tmin'] = _np.min(self.data[:,:,:,:,3])
        self.header['tmax'] = _np.max(self.data[:,:,:,:,3])
        self.header['nt']   = _np.shape(self.data)[0]
