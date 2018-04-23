import tarfile 
from Tfs import Cast 
import numpy as _numpy
import copy as _copy

class TfsArray : 

    def __init__(self,filename=None) : 
        self.header  = {}
        self.columns = [] 
        self.formats = [] 
        self.segment = [] 
        self.data    = []
        self.nitems  = 0
        self.filename= filename
        if self.filename != None : 
            self.Load(self.filename) 


    def Clear(self) :
        """ 
        Clear() 
        empty all internal variables 
        """
        self.__init__()

    def Load(self,filename) :
        """
        Load('filename.tfs')        
        read the tfs file and prepare data structures
        """

        # open gz file 
        if 'gz' in filename :
            print 'pymadx.TfsArray.Load> zipped file'
            tar = tarfile.open(filename,'r')
            f = tar.extractfile(tar.firstmember)
        else:
            print 'pymadx.TfsArray> normal file'
            f = open(filename)


        self.isegment = -1
        self.columns.append('SEGMENT')
        self.formats.append('%d')

        for line in f :
            sl = line.strip('\n').split()
            if line[0] == "@" : 
                # header
                # print 'pymadx.TfsArray> header',line 
                self.header[sl[1]] = sl[-1].replace('"','')
            elif line[0] == '*' : 
                # column name 
                # print 'pymadx.TfsArray> columns',line 
                self.columns.extend(sl[1:]) #miss *                
            elif line[0] == '$' : 
                # column type
                # print 'pymadx.TfsArray> type',line 
                self.formats.extend(sl[1:]) #miss $              
            elif sl[0] == '#segment' : 
                # entry marker 
                # print 'pymadx.TfsArray> marker',line, self.nitems,isegment
                self.isegment=self.isegment+1
                
            else :                 
                # data
                d = [Cast(item) for item in sl]
                d.insert(0,self.isegment) 
                self.data.append(d)    # put in data dict by name
                self.nitems += 1 # keep tally of number of items                

            
        # make numpy array for convenience
        self.dataArray = _numpy.array(self.data)

    def __repr__(self):
        s =  ''
        s += 'pymadx.TfsArray instance\n'
        s += str(self.nitems) + ' items in lattice\n'
        return s

    def __iter__(self):
        self._iterindex = -1
        return self

    def next(self):
        if self._iterindex == len(self.sequence)-1:
            raise StopIteration
        self._iterindex += 1
        return self.GetRow(self._iterindex)

    def __getitem__(self,index):
        #return single item or slice of lattice
        if type(index) == slice:
            #prepare step integer - allow reverse stepping too
            if index.stop > index.start:
                index = slice(index.start,index.stop,1)
            else:
                index = slice(index.start,index.stop,-1)
            return [self.GetRow(i) for i in range(index.start,index.stop,index.step)]
        elif type(index) == int:
            return self.GetRow(index)
        else:
            raise ValueError("argument not an index or a slice")

    def ColumnIndex(self,columnstring):
        """
        ColumnIndex(columnname):
        
        return the index to the column matching the name
        
        REMEMBER: excludes the first column NAME
        0 counting

        """
        return self.columns.index(columnstring)

    def NameFromIndex(self,index):
        return self.dataArray[index,self.columns.index['NAME']]
            
    def GetColumn(self, colname) :
        colindex = self.columns.index(colname) 
        return self.dataArray[:,colindex]

    def GetRow(self, rowindex) : 
        return self.dataArray[rowindex,:]

    def GetSegment(self, segmentindex) : 
        segcut = self.dataArray[:,0] == segmentindex
        c = _copy.deepcopy(self)
        c.data      = self.data
        c.dataArray = self.dataArray[segcut,:]
        c.nitems    = len(c.dataArray)        
        return c

    def InterrogateItem(self,index):
        """
        InterrogateItem(index)
        
        Print out all the parameters and their
        names for a particlular element in the 
        sequence identified by name
        """
        for i,parameter in enumerate(self.columns):
            print parameter.ljust(10,'.'),self.dataArray[index,i]

    def GetElementsOfType(self,typename) :
        """
        GetElementsOfType(typename) 
        
        Returns a list of the names of elements of a certain type

        typename can be a single string or a tuple or list of strings

        ie 
        GetElementsOfType('SBEND')
        GetElementsOfType(['SBEND','RBEND'])
        GetElementsOfType(('SBEND','RBEND','QUADRUPOLE'))

        """        
        i = self.ColumnIndex('KEYWORD')
        types = self.GetColumn('KEYWORD')
        return [t for t in types if t in typename]

    def ReportPopulations(self):
        """
        Print out all the population of each type of
        element in the beam line (sequence)
        """
        print 'Filename > ',self.filename
        print 'Total number of items > ',self.nitems
        i = self.ColumnIndex('KEYWORD')
        types = set(self.GetColumn('KEYWORD'))
        populations = [(len(self.GetElementsOfType(key)),key) for key in types]
        print 'Type'.ljust(15,'.'),'Population'
        for item in sorted(populations)[::-1]:
            print item[1].ljust(15,'.'),item[0]
