# pybdsim.Builder - tools to build bdsim lattices
# Version 1.0
# L. Nevay
# laurie.nevay@rhul.ac.uk

"""
Writer

Write files for a pybdsim.Builder.Machine instance. Each section of the written output
(e.g. components, sequence, beam etc.) can be written in the main gmad file, written 
in its own separate file, or called from an external, pre-existing file.

Classes:
File - A class that represents each section of the written output - contains booleans and strings.
Writer - A class that writes the data to disk.

"""
from . import _General
from . import Beam as _Beam
from . import Options as _Options
from . import Builder as _Builder
import time as _time
import os as _os
import numpy as _np
import textwrap as _textwrap
sections = ['components',
            'sequence',
            'samplers',
            'beam',
            'options',
            'bias']

class FileSection():
    """
    A class that represents a section of a gmad file. The sections that this
    class can represent are:
    
    * Components
    * Sequence
    * Samplers
    * Beam
    * Options
    * Bias

    The class contains booleans and strings relating to the location of 
    that sections data. The section can set to be:
    
    * Written in its own separate file (default)
    * Written in the main gmad file
    * Called from an external file

    These classes are instantiated in the writer class for each section. 
    An optional string passed in upon class instantiation is purely for 
    the representation of the object which will state where the data will 
    be written/called. This string should be one of the section names 
    listed above.
    
    Example:

    >>> beam = FileSection('beam')
    >>> beam.CallExternalFile('../myBeam.gmad')
    >>> beam
    pybdsim.Writer.File instance
    File data will be called from the external file:
    ../myBeam.gmad  
    
    """
    
    def __init__(self,willContain=''):
        #bool for file location/calling
        self._writeInMain         = False
        self._isWrittenSeparately = True
        self._isUserDefined       = False
        self._filePath            = ''
        
        #bool for written status, used to determine if required sections have been written
        self._hasBeenWritten      = False
        
        #string for __repr__ output
        self._willContain=''
        if isinstance(willContain, str) and willContain in sections:
            self._willContain = willContain

    def __repr__(self):
        s = ''
        s += 'pybdsim.Writer.File instance\n'
        
        if self._writeInMain:
            s += 'File data will be written into the main file.\n'
        elif self._isWrittenSeparately:
            if self._willContain != '':
                self._willContain += ' '
            s += 'File data will be written into a separate '+self._willContain+'file.\n'
        elif self._isUserDefined:
            s += 'File data will be called from the external file:\n'
            s += (self._filePath+'\n')
        return s

    def CallExternalFile(self,filepath=''):
        if not isinstance(filepath,_np.str):
            raise TypeError("Filepath must be a string.")
        self._isUserDefined       = True
        self._isWrittenSeparately = False
        self._writeInMain         = False
        self._filePath            = filepath
        self._hasBeenWritten      = True

    def WriteInMain(self):
        self._isUserDefined       = False
        self._isWrittenSeparately = False
        self._writeInMain         = True
        self._filePath            = ''

    def WriteSeparately(self):
        self._isUserDefined       = False
        self._isWrittenSeparately = True
        self._writeInMain         = False
        self._filePath            = ''


class Writer():
    """
    A class for writing a pybdsim.Builder.Machine instance to file.
    
    This class allows the user to write individual sections of a BDSIM input file
    (e.g. components, sequence, beam etc.) or write the machine as a whole.
    
    There are 6 attributes in this class which are FileSection instances representing each
    section of the data. The location where these sections will be written/read is
    stored in these instances. See the FileSection class for further details. 
    
    The optional boolean 'singlefile' in the WriteMachine function for writing
    the sections to a single file overrides any sections locations set in their
    respective FileSection instances.
    
    This class also has individual functions (e.g. WriteBeam) to write each file section
    and the main file (WriteMain) separately. These section functions must be called BEFORE
    the WriteMain function is called otherwise the main file will have no reference to 
    these sections.
    
    Examples:
    
    Writing the Builder.Machine instance myMachine to separate files:
    
    >>> a = Writer()
    >>> a.WriteMachine(myMachine,'lattice.gmad')
    Lattice written to:
    lattice_components.gmad
    lattice_sequence.gmad
    lattice_beam.gmad
    lattice.gmad
    All included in main file:
    lattice.gmad

    Writing the Builder.Machine instance myMachine into a single file:

    >>> a = Writer()
    >>> a.WriteMachine(myMachine,'lattice.gmad',singlefile=True)
    Lattice written to:
    lattice.gmad
    All included in main file:
    lattice.gmad

    """
    def __init__(self):
        #FileSection instances for each section of the input.
        self.Components = FileSection('components')
        self.Sequence   = FileSection('sequence')
        self.Samplers   = FileSection('samplers')
        self.Beam       = FileSection('beam')
        self.Options    = FileSection('options')
        self.Bias       = FileSection('bias')
        
        self._mainFileLines   = []        #lines that will be written to the main file.
        self._elementsperline = 100       #number of machine elements per bdsim line (not text line)
        self._basefilename    = 'lattice' #default base file name
        self._timestring      = '! ' + _time.strftime("%a, %d %b %Y %H:%M:%S +0000", _time.gmtime()) + '\n'

        #default filenames. These can be overwritten to a user defined filename.
        self._mainFilename = self._basefilename + '.gmad'
        self._defaultSectionFilenames  = {}
        for sectiontype in sections:
            self._defaultSectionFilenames[sectiontype] = self._basefilename + '_' + sectiontype + '.gmad'
        
        #list of sections that will be / have been written. This is for both user feedback and
        #for any include lines that will be written in the main file.
        self._sectionsToBeWritten = []

    def WriteMachine(self,machine, filename, singlefile=False,
                     verbose=True, summary=True, overwrite=True):
        """
        WriteMachine(machine(machine),filename(string),singlefile(bool),verbose(bool))
        
        Write a machine to disk. By default, the machine will be written
        into the following individual files:
        
        +---------------------------+----------------------------------------+
        | filename_components.gmad  | component files (max 10k per file)     |
        +---------------------------+----------------------------------------+
        | filename_sequence.gmad    | lattice definition                     |
        +---------------------------+----------------------------------------+
        | filename_samplers.gmad    | sampler definitions (max 10k per file) |
        +---------------------------+----------------------------------------+
        | filename_options.gmad     | options                                |
        +---------------------------+----------------------------------------+
        | filename_beam.gmad        | beam definition                        |
        +---------------------------+----------------------------------------+
        | filename_bias.gmad        | machine biases (if defined)            |
        +---------------------------+----------------------------------------+
        | filename.gmad             | suitable main file with all sub        |
        |                           | files in correct order                 |
        +---------------------------+----------------------------------------+

        These are prefixed with the specified filename / path

        The optional bool singlefile = True will write all the above sections
        into a single file:

        filename.gmad

        kwargs:
        overwrite : Do not append an integer to the basefilename if
        already exists, instead overwrite existing files.

        """
        self._checkFiles(filename, overwrite)

        if singlefile:
            #set all sections to write in the main file
            self.Components.WriteInMain()
            self.Sequence.WriteInMain()
            self.Beam.WriteInMain()
            self.Options.WriteInMain()
            self.Samplers.WriteInMain()
            self.Bias.WriteInMain()

        #write the individual files.
        self.WriteComponents(machine)
        self.WriteSequence(machine)
        self.WriteSamplers(machine)
        self.WriteBeam(machine)
        self.WriteOptions(machine)
        self.WriteBias(machine)
        
        #Write main
        self.WriteMain(machine, summary=summary)

        if verbose:
            #user feedback
            print('Lattice written to:')
            for section in self._sectionsToBeWritten:
                sectObject = getattr(self,section) #a copy of the FileSection object.
                fn = getattr(sectObject,'_filePath')
                print(fn)
            print('All included in main file: \n',self._mainFilename)

    def WriteMain(self,machine,filename='',summary=True):
        """
        WriteMain(machine(machine),filename(string))
        
        Write the main gmad file:
        filename.gmad
        
        The functions for the other sections of the machine (components,sequence,beam,options,samplers,bias)
        must be written BEFORE this function is called.
        
        """
        self._machineCheck(machine)
        
        if filename == '':
            fn_main = self._mainFilename #default
        else:
            fn_main = filename #override

        #do not write if the beam, components or sequence sections have not been written.
        #These three have to exist for a machine to work.
        if (not self.Components._hasBeenWritten) or (not self.Sequence._hasBeenWritten) or (not self.Beam._hasBeenWritten):
            exceptionString = "The following sections must be written before the main file can be written:\n"
            exceptionString += "- Components\n- Sequence\n- Beam"
            raise Exception(exceptionString)

        #write main file
        f = open(fn_main,'w')
        f.write(self._timestring)
        if (summary):
            f.write('! pybdsim.Builder Lattice \n')
            f.write('! number of elements = ' + str(len(machine.elements)) + '\n')
            f.write('! total length       = ' + str(machine.length) + ' m\n\n')
        else:
            f.write('\n\n')
            
        #other files to include
        for section in self._sectionsToBeWritten:
            sectObject = getattr(self,section)
            fn = getattr(sectObject,'_filePath')
            isUserDefined = getattr(sectObject,'_isUserDefined')
            if not isUserDefined:
                fn = fn.split('/')[-1]
            f.write('include '+fn+';\n')
        f.write('\n\n')
        
        # write lines to main file from components, sequence etc.
        for line in self._mainFileLines:
            f.write(line)
        #write samplers in this main file if less than 10 samplers.
        if len(machine.samplers) <= 10:
            for sampler in machine.samplers:
                f.write(str(sampler))
        f.close()
        self._mainFilename = fn_main

    def WriteComponents(self,machine,filename=''):
        """
        Write the machines components to disk:
        filename.gmad
        """
        self._machineCheck(machine)
        fn_components = self._getName(filename,'components')
        if self.Components._writeInMain:                #if _writeInMain, append strings to _mainFileLines list.
            for element in machine.elements:
                self._mainFileLines.append(str(element))
            self._mainFileLines.append('\r\n')
        elif self.Components._isWrittenSeparately:      #if _isWrittenSeparately, write directly to file here.
            f = open(fn_components,'w')
            f.write(self._timestring)
            f.write('! pybdsim.Builder Lattice \n')
            f.write('! COMPONENT DEFINITION\n\n')
            for element in machine.elements:
                f.write(str(element))
            f.close()
            self._sectionsToBeWritten.append('Components')
            self.Components._filePath = fn_components   #update FileSection path
        elif self.Components._isUserDefined:
            self._sectionsToBeWritten.append('Components')
        self.Components._hasBeenWritten = True

    def WriteBias(self,machine,filename=''):
        """
        Write the machines bias to disk:
        filename.gmad
        """
        self._machineCheck(machine)
        fn_bias = self._getName(filename,'bias')

        #write bias if it exists
        if len(machine.bias) > 0:
            if self.Bias._writeInMain:
                for bias in machine.bias:
                    self._mainFileLines.append(str(bias))
                self._mainFileLines.append('\r\n')
            elif self.Bias._isWrittenSeparately:
                f = open(fn_bias,'w')
                for bias in machine.bias:
                    f.write(str(bias))
                f.close()
                self.Bias._filePath = fn_bias
                self._sectionsToBeWritten.append('Bias')
        elif self.Bias._isUserDefined:
            self._sectionsToBeWritten.append('Bias')
        self.Bias._hasBeenWritten = True

    def WriteBeam(self,machine,filename=''):
        """
        Write a machines beam to disk:
        filename.gmad
        
        Machine can be either a pybdsim.Builder.Machine instance
        or a pybdsim.Beam.Beam instance.
        """
        #check data type
        if isinstance(machine,_Beam.Beam):
            object = machine
        else:
            self._machineCheck(machine)
            object = machine.beam
        
        fn_beam = self._getName(filename,'beam')

        if self.Beam._writeInMain:
            self._mainFileLines.append(object.ReturnBeamString())
            self._mainFileLines.append('\r\n')
        elif self.Beam._isWrittenSeparately:
            f = open(fn_beam,'w')
            f.write(self._timestring)
            f.write('! pybdsim.Builder \n')
            f.write('! BEAM DEFINITION \n\n')
            f.write(object.ReturnBeamString())
            f.close()
            self.Beam._filePath = fn_beam
            self._sectionsToBeWritten.append('Beam')
        elif self.Beam._isUserDefined:
            self._sectionsToBeWritten.append('Beam')
        self.Beam._hasBeenWritten = True

    def WriteSamplers(self,machine,filename=''):
        """
        Write the machines samplers to disk:
        filename.gmad
        """
        self._machineCheck(machine)
        fn_samplers = self._getName(filename,'samplers')

        #write samplers
        # if less than 10 samplers, just put in main file
        if len(machine.samplers) > 10:
            if self.Samplers._writeInMain:
                for sampler in machine.samplers:
                    self._mainFileLines.append(str(sampler))
                self._mainFileLines.append('\r\n')
            elif self.Samplers._isWrittenSeparately:
                f = open(fn_samplers,'w')
                f.write(self._timestring)
                f.write('! pybdsim.Builder \n')
                f.write('! SAMPLER DEFINITION\n\n')
                for sampler in machine.samplers:
                    f.write(str(sampler))
                f.close()
                self.Samplers._filePath = fn_samplers
                self._sectionsToBeWritten.append('Samplers')
        elif self.Samplers._isUserDefined:
            self._sectionsToBeWritten.append('Samplers')
        self.Samplers._hasBeenWritten = True

    def WriteOptions(self,machine,filename=''):
        """
        Write a machines options to disk:
        filename.gmad
        
        Machine can be either a pybdsim.Builder.Machine instance
        or a pybdsim.Options.Options instance.
        """
        #check data type
        if isinstance(machine,_Options.Options):
            object = machine
        else:
            self._machineCheck(machine)
            object = machine.options
        
        fn_options = self._getName(filename,'options')

        # write options - only if specified
        if object != None:
            if self.Options._writeInMain:
                self._mainFileLines.append(object.ReturnOptionsString())
                self._mainFileLines.append('\r\n')
            elif self.Options._isWrittenSeparately:
                f = open(fn_options,'w')
                f.write(self._timestring)
                f.write('! pybdsim.Builder \n')
                f.write('! OPTIONS DEFINITION \n\n')
                f.write(object.ReturnOptionsString())
                f.close()
                self.Options._filePath = fn_options
                self._sectionsToBeWritten.append('Options')
        elif self.Options._isUserDefined:
            self._sectionsToBeWritten.append('Options')
        self.Options._hasBeenWritten = True

    def WriteSequence(self,machine,filename=''):
        """
        Write the machines sequence to disk:
        filename.gmad
        """
        self._machineCheck(machine)
        fn_sequence = self._getName(filename,'sequence')

        # need to define the period before making sampler planes
        sequencechunks = _General.Chunks(machine.sequence,self._elementsperline)
        linelist = []
        ti = 0

        #write lattice sequence
        if self.Sequence._writeInMain:
            for line in _General.Chunks(machine.sequence,self._elementsperline):
                # Use _textwrap.wrap to wrap very long lines
                linetxt = '\n\t'.join(_textwrap.wrap(
                    "l{}: line = ({});".format(ti, ', '.join(line))))
                self._mainFileLines.append("{}\n".format(linetxt))
                linelist.append('l'+str(ti))
                ti += 1
            self._mainFileLines.append('lattice: line = ('+', '.join(linelist)+');\n')
            self._mainFileLines.append('use, period=lattice;\n')
            self._mainFileLines.append('\r\n')
        elif self.Sequence._isWrittenSeparately:
            f = open(fn_sequence,'w')
            f.write(self._timestring)
            f.write('! pybdsim.Builder \n')
            f.write('! LATTICE SEQUENCE DEFINITION\n\n')
            for line in _General.Chunks(machine.sequence,self._elementsperline):
                # Use _textwrap.wrap to wrap very long lines
                linetxt = '\n\t'.join(_textwrap.wrap(
                    "l{}: line = ({});".format(ti, ', '.join(line))))
                f.write("{}\n".format(linetxt))
                linelist.append('l'+str(ti))
                ti += 1

            f.write('lattice: line = ('+', '.join(linelist)+');\n')
            f.write('use, period=lattice;\n')
            f.close()
            self.Sequence._filePath = fn_sequence
            self._sectionsToBeWritten.append('Sequence')
        elif self.Sequence._isUserDefined:
            self._sectionsToBeWritten.append('Sequence')
        self.Sequence._hasBeenWritten = True

    def _getName(self,filename,sectiontype=''):
        #check input types are strings
        if not isinstance(filename,_np.str):
            raise TypeError("Filename not a string")
        if not isinstance(sectiontype,_np.str):
            raise TypeError("Sectiontype not a string")
        
        if filename == '' and sectiontype == '':
            raise ValueError("Both filename and sectiontype cannot be empty strings")
        
        #get section file name
        if sectiontype in sections:
            if (filename == '') or (filename == self._basefilename):
                fn_name = self._defaultSectionFilenames[sectiontype]    #no name or name = basefilename, use section default
            else:
                fn_name = self._checkExtensionAndPath(filename)         #otherwise, override the filename
        else:
            if (filename == '') or (filename == self._basefilename):    #no section but filename clashes with basefilename
                errorString =  "Unknown section type - filename cannot be an empty string "
                errorString += "or match the basefilename ("+_np.str(self._basefilename)+")."
                raise ValueError(errorString)
            else:
                fn_name = self._checkExtensionAndPath(filename)         #override filename

        #check for duplicate file names.
        for section in self._sectionsToBeWritten:
            sectObject = getattr(self,section)
            fn = getattr(sectObject,'_filePath')
            if (fn == fn_name) or (fn == self._mainFilename):
                raise ValueError("Filename already used.")
        return fn_name

    def _machineCheck(self,machine):
        if not isinstance(machine,_Builder.Machine):
            raise TypeError("Not a machine instance")

    def _checkFiles(self,filename, overwrite=True):
        filename = self._checkExtensionAndPath(filename)

        if not overwrite:
            #check if file already exists
            originalFilename = filename
            filename = _General.GenUniqueFilename(filename)
            if filename != originalFilename:
                print('Warning, chosen filename already exists - using filename: ',filename.split('.')[0])

        basefilename = filename[:-5] #everything before '.gmad'
        #new default section names
        self._defaultSectionFilenames['components'] = basefilename + '_components.gmad'
        self._defaultSectionFilenames['sequence']   = basefilename + '_sequence.gmad'
        self._defaultSectionFilenames['samplers']   = basefilename + '_samplers.gmad'
        self._defaultSectionFilenames['beam']       = basefilename + '_beam.gmad'
        self._defaultSectionFilenames['options']    = basefilename + '_options.gmad'
        self._defaultSectionFilenames['bias']       = basefilename + '_bias.gmad'
        self._mainFilename = basefilename + '.gmad'
        self._basefilename = basefilename
        self._timestring = '! ' + _time.strftime("%a, %d %b %Y %H:%M:%S +0000", _time.gmtime()) + '\n'

    def _checkExtensionAndPath(self,filename):
        #check filename for extension
        if filename[-5:] != '.gmad':
            filename += '.gmad'

        #check for directory and make it if not:
        if '/' in filename:
            directory = '/'.join(filename.split('/')[:-1]) #strip the filename off
            if not _os.path.exists(directory):
                _os.system("mkdir -p " + directory)
        return filename

