import numpy as _np

from pybdsim import Builder as _pyBuilder
from pybdsim import Options as _Options

from pytransport.Elements import Elements as _Elements
from pytransport import Reader as _Reader
from pytransport import _General
from pytransport.Data import ConversionData


def Transport2Gmad(inputfile,
                   particle='proton',
                   distrType='gauss',
                   outputDir='gmad',
                   debug=False,
                   dontSplit=False,
                   keepName=False,
                   combineDrifts=False):
    """
    **Transport2Gmad** convert a Transport input or output file into a gmad input file for bdsim

    +-------------------------------+-------------------------------------------------------------------+
    | **inputfile**                 | dtype = string                                                    |
    |                               | path to the input file                                            |
    +-------------------------------+-------------------------------------------------------------------+
    | **particle**                  | dtype = string. Optional, default = "proton"                      |
    |                               | the particle species                                              |
    +-------------------------------+-------------------------------------------------------------------+
    | **distrType**                 | dtype = string. Optional, Default = "gauss".                      |
    |                               | the beam distribution type. Can be either gauss or gausstwiss.    |
    +-------------------------------+-------------------------------------------------------------------+
    | **outputDir**                 | dtype=string. Optional, default = "gmad"                          |
    |                               | the output directory where the files will be written              |
    +-------------------------------+-------------------------------------------------------------------+
    | **debug**                     | dtype = bool. Optional, default = False                           |
    |                               | output a log file (inputfile_conversion.log) detailing the        |
    |                               | conversion process, element by element                            |
    +-------------------------------+-------------------------------------------------------------------+
    | **dontSplit**                 | dtype = bool. Optional, default = False                           |
    |                               | the converter splits the machine into multiple parts when a beam  |
    |                               | is redefined in a Transport lattice. dontSplit overrides this and |
    |                               | forces the machine to be written to a single file                 |
    +-------------------------------+-------------------------------------------------------------------+
    | **keepName**                  | dtype = bool. Optional, default = False                           |
    |                               | keep the names of elements as defined in the Transport inputfile  |
    +-------------------------------+-------------------------------------------------------------------+
    | **combineDrifts**             | dtype = bool. Optional, default = False                           |
    |                               | combine multiple consecutive drifts into a single drift           |
    +-------------------------------+-------------------------------------------------------------------+

    Example:

    >>> Transport2Gmad(inputfile)

    Writes converted machine to disk. Reader automatically detects if the supplied input file is a Transport input
    file or Transport output file.

    """
    # Explicitly specify. Auto is set to True to automatically convert, but is not an input argument. No madx conversion.
    _Convert(inputfile=inputfile, particle=particle, debug=debug, distrType=distrType, gmadDir=outputDir,
             auto=True, dontSplit=dontSplit, keepName=keepName, combineDrifts=combineDrifts, outlog=True)


class _Convert(_Elements):
    def __init__(self, inputfile, particle, debug, distrType, gmadDir, auto, dontSplit, keepName, combineDrifts, outlog):
        # instantiate the main data container and pass in to element class.
        # Explicitly specify as correct machine needs passing in.
        transportData = ConversionData(inputfile=inputfile, options=_Options.Options(), machine=_pyBuilder.Machine(),
                                       particle=particle, debug=debug, distrType=distrType, gmad=True, gmadDir=gmadDir,
                                       madx=False, madxDir='', auto=auto, dontSplit=dontSplit, keepName=keepName,
                                       combineDrifts=combineDrifts, outlog=outlog)
        _Elements.__init__(self, transportData)
        # convert automatically
        if auto:
            self.Transport2Gmad()

    def LoadFile(self, inputfile):
        """
        Load a Transport file.
        """
        temp = _Reader.Reader()

        if not isinstance(inputfile, _np.str):
            raise TypeError("Input must be a string")

        infile = inputfile.split('/')[-1]  # Remove filepath, leave just filename
        self.Transport._file = infile[:-4]  # Remove extension
        self.Transport._filename = inputfile
        isOutput = _General.CheckIsOutput(inputfile)  # Is a TRANSPORT standard output file.
        self.Writer.DebugPrintout("File Read.")
        if isOutput:
            lattice, output = temp.GetLatticeAndOptics(inputfile)
            fits, fitres = temp.GetFits(inputfile)
            self.Transport = _General.OutputFitsToRegistry(self.Transport, fitres)
            self.Writer.DebugPrintout('Adding any fitting output to the fitting registry (self.FitRegistry)')

            self.Writer.DebugPrintout('Processing file and adding to Transport class.')
            for linenum, latticeline in enumerate(lattice):
                latticeline = latticeline.replace(';', '')
                line = _np.array(latticeline.split(' '), dtype=_np.str)
                line = _General.RemoveIllegals(line)

                # Method of dealing with split lines in the output
                # Should only be applicable to type 12 entry (up to 15 variables)
                # It is assumed that the line is always split, so be careful.
                prevline = lattice[linenum - 1].replace(';', '')
                prevline = _np.array(prevline.split(' '), dtype=_np.str)
                prevline = _General.RemoveIllegals(prevline)

                try:
                    if (linenum > 0) and _np.abs(_np.float(line[0])) == 12.0:
                        latticeline, line = _General.JoinSplitLines(linenum, lattice)
                    # Ignore line after type 12 entry (second part of split line)
                    if (linenum > 1) and _np.abs(_np.float(prevline[0])) == 12.0:
                        pass
                    else:
                        self.Transport.data.append(line)
                        self.Transport.filedata.append(latticeline)
                except ValueError:
                    self.Transport.data.append(line)
                    self.Transport.filedata.append(latticeline)
                except IndexError:
                    pass

        else:
            f = open(inputfile)
            for inputline in f:
                endoflinepos = _General.FindEndOfLine(inputline)
                templine = inputline
                if endoflinepos > 0:
                    templine = inputline[:endoflinepos]
                line = _np.array(templine.split(' '), dtype=_np.str)
                # do not change comment lines
                if not line[0][0] == '(':
                    line = _General.RemoveIllegals(line)
                self.Transport.data.append(line)
                self.Transport.filedata.append(inputline)
            f.close()
        self.Transport.convprops.fileloaded = True

    def Write(self):
        """
        Write the converted TRANSPORT file to disk.
        """
        self.Writer.DebugPrintout("Adding beam to gmad machine:")
        self.Transport.AddBeam()
        self.Transport.AddOptions()
        self.Transport.machine.AddSampler('all')
        self.Writer.BeamDebugPrintout(self.Transport.beamprops, self.Transport.units)
        fname = _General.RemoveFileExt(self.Transport.convprops.file)
        if self.Transport.convprops.numberparts < 0:
            filename = fname
        else:
            self.Transport.convprops.numberparts += 1
            filename = fname + '_part' + _np.str(self.Transport.convprops.numberparts)
        self.Writer.Write(self.Transport, filename)

    def Transport2Gmad(self):
        """
        Convert, process, and write.
        """
        self.LoadFile(self.Transport.convprops.file)

        if not self.Transport.convprops.fileloaded:
            self.Writer.Printout('No file loaded.')
            return
        self.ProcessAndBuild()
        self.Write()

    def _ElementPrepper(self, line, linenum, filetype='input'):
        """
        Function to extract the data and prepare it for processing by each element function.
        This has been written as the lattice lines from an input file and output file are different,
        so it just a way of correctly ordering the information.
        """
        linedict = {'elementnum': 0.0,
                    'name': '',
                    'length': 0.0,
                    'isZeroLength': True}
        numElements = _np.str(len(self.Transport.ElementRegistry.elements))
        typeNum = _General.GetTypeNum(line)
        linedict['elementnum'] = typeNum

        if typeNum == 15.0:
            label = _General.GetLabel(line)
            if filetype == 'output':
                linedict['label'] = line[2].strip('"')
            if filetype == 'input':
                linedict['label'] = label
            linedict['number'] = line[1]
            self.Writer.ElementPrepDebugPrintout("Unit Control", numElements)

        if typeNum == 20.0:
            angle = 0  # Default
            if len(line) == 2:  # i.e has a label
                endofline = _General.FindEndOfLine(line[1])
                angle = line[1][:endofline]
            else:
                for index in line[1:]:
                    try:
                        angle = _np.float(index)
                        break
                    except ValueError:
                        pass
            linedict['angle'] = angle
            self.Writer.ElementPrepDebugPrintout("coordinate rotation", numElements)

        if typeNum == 1.0:
            linedict['name'] = _General.GetLabel(line)
            linedict['isAddition'] = False
            if _General.CheckIsAddition(line, filetype):
                linedict['isAddition'] = True
            # line = self._remove_label(line)
            if len(line) < 8:
                raise IndexError("Incorrect number of beam parameters.")
            n = 1
            if filetype == 'input':
                n = 0
            elif filetype == 'output':
                n = 1

            # Find momentum
            # endofline = self._endofline(line[7+n])
            # if endofline != -1:
            #     linedict['momentum'] = line[7+n][:endofline]
            # else:
            linedict['momentum'] = line[7 + n]
            linedict['Sigmax'] = line[1 + n]
            linedict['Sigmay'] = line[3 + n]
            linedict['Sigmaxp'] = line[2 + n]
            linedict['Sigmayp'] = line[4 + n]
            linedict['SigmaT'] = line[5 + n]
            linedict['SigmaE'] = line[6 + n]
            self.Writer.ElementPrepDebugPrintout("Beam definition or r.m.s addition", numElements)

        if typeNum == 2.0:
            linedict['name'] = _General.GetLabel(line)
            linedict['data'] = _General.GetElementData(line)
            self.Writer.ElementPrepDebugPrintout("poleface rotation", numElements)

        if typeNum == 3.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("drift", numElements)

        if typeNum == 4.0:
            linedict['name'] = _General.GetLabel(line)
            linedict['linenum'] = linenum
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            e1, e2 = _General.GetFaceRotationAngles(self.Transport.data, linenum)
            linedict['e1'] = e1
            linedict['e2'] = e2
            self.Writer.ElementPrepDebugPrintout("dipole", numElements)

        if typeNum == 5.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("quadrupole", numElements)

        if typeNum == 6.0:
            # element is a collimator or transform update
            # transform update is later ignored so only update linedict as if collimator

            physicalElements = [1.0, 3.0, 4.0, 5.0, 11.0, 18.0, 19.0]

            # Only iterate if not the last element
            if linenum == len(self.Transport.data):
                pass
            else:
                # Since collimators have zero length in TRANSPORT, chosen to use length of next drift instead if
                # present. Check all remaining elements for the next drift, following element(s) may be non-physical
                # element which can be ignored as it shouldnt affect the beamline. Ignore beam definition too in
                # the case where machine splitting is not permitted.
                for nextline in self.Transport.data[linenum + 1:]:
                    nextTypeNum = _General.GetTypeNum(nextline)
                    if nextTypeNum == 3.0:
                        nextData = _General.GetElementData(nextline)
                        linedict['length'] = nextData[0]
                        linedict['isZeroLength'] = False
                        linedict['name'] = _General.GetLabel(line)
                        data = _General.GetElementData(line)
                        linedict['data'] = data
                        break
                    # stop if physical element or beam redef if splitting permitted
                    elif nextTypeNum in physicalElements:
                        if (nextTypeNum == 1.0) and self.Transport.convprops.dontSplit:
                            pass
                        elif (nextTypeNum == 6.0) and self.Transport.convprops.typeCode6IsTransUpdate:
                            pass
                        else:
                            break
                    # ignore non-physical element
                    else:
                        pass
            # Can be either transform update or collimator, a 16. 14 entry changes the definition but is only
            # processed after ALL elements are added to the registry.
            self.Writer.ElementPrepDebugPrintout("transform update or collimator", numElements)

        if typeNum == 9.0:
            self.Writer.ElementPrepDebugPrintout("repetition control", numElements)

        if typeNum == 11.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['voltage'] = data[1]
            linedict['isZeroLength'] = False
            if len(data) == 4:  # Older case for single element
                linedict['phase_lag'] = data[2]
                linedict['wavel'] = data[3]
            self.Writer.ElementPrepDebugPrintout("acceleration element", numElements)

        if typeNum == 12.0:
            linedict['data'] = _General.GetElementData(line)
            linedict['name'] = _General.GetLabel(line)

            prevline = self.Transport.data[linenum - 1]  # .split(' ')
            linedict['prevlinenum'] = _np.float(prevline[0])
            linedict['isAddition'] = False
            if _General.CheckIsAddition(line):
                linedict['isAddition'] = True
            self.Writer.ElementPrepDebugPrintout("beam rotation", numElements)

        if typeNum == 13.0:
            linedict['data'] = _General.GetElementData(line)
            self.Writer.ElementPrepDebugPrintout("Input/Output control", numElements)

        if typeNum == 16.0:
            linedict['data'] = _General.GetElementData(line)
            self.Writer.ElementPrepDebugPrintout("special input", numElements)

        if typeNum == 18.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("sextupole", numElements)

        if typeNum == 19.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("solenoid", numElements)

        if typeNum == 22.0:
            self.Writer.ElementPrepDebugPrintout("space charge element", numElements)

        if typeNum == 23.0:
            self.Writer.ElementPrepDebugPrintout("buncher", numElements)

        rawline = self.Transport.filedata[linenum]
        self.Transport.ElementRegistry.AddToRegistry(linedict, rawline)

    def ProcessAndBuild(self):
        """
        Function that loops over the lattice, adds the elements to the element registry,
        and updates any elements that have fitted parameters.
        It then converts the registry elements and adds to the gmad machine.
        """
        self.Writer.DebugPrintout('Processing tokenised lines from input file and adding to element registry.\n')

        for linenum, line in enumerate(self.Transport.data):
            self.Writer.DebugPrintout('Processing tokenised line ' + _np.str(linenum) + ' :')
            self.Writer.DebugPrintout('\t' + str(line))
            self.Writer.DebugPrintout('\tOriginal :')
            self.Writer.DebugPrintout('\t' + self.Transport.filedata[linenum])

            # Checks if the SENTINEL line is found. SENTINEL relates to TRANSPORT fitting routine and is only written
            # after the lattice definition, so there's no point reading lines beyond it.
            if _General.CheckIsSentinel(line):
                self.Writer.DebugPrintout('Sentinel Found.')
                break
            # Test for positive element, negative ones ignored in TRANSPORT so ignored here too.
            try:
                typeNum = _General.GetTypeNum(line)
                if typeNum > 0:
                    if self.Transport.data[0][0] == 'OUTPUT':
                        self._ElementPrepper(line, linenum, 'output')
                        self._UpdateElementsFromFits()
                    else:
                        line = _General.RemoveIllegals(line)
                        self._ElementPrepper(line, linenum, 'input')
                else:
                    self.Writer.DebugPrintout('\tType code is 0 or negative, ignoring line.')
            except ValueError:
                errorline = '\tCannot process line ' + _np.str(linenum) + ', '
                if line[0][0] == '(' or line[0][0] == '/':
                    errorline += 'line is a comment.'
                elif line[0][0] == 'S':  # S used as first character in SENTINEL command.
                    errorline = 'line is for TRANSPORT fitting routine.'
                elif line[0] == '\n':
                    errorline = 'line is blank.'
                else:
                    errorline = 'reason unknown.'
                self.Writer.DebugPrintout(errorline)

        self.Writer.DebugPrintout(
            'Converting registry elements to pybdsim compatible format and adding to machine builder.\n')

        skipNextDrift = False  # used for collimators
        lastElementWasADrift = True  # default value
        if self.Transport.convprops.combineDrifts:
            lastElementWasADrift = False
        for linenum, linedict in enumerate(self.Transport.ElementRegistry.elements):
            self.Writer.DebugPrintout('Converting element number ' + _np.str(linenum) + ':')
            convertline = '\t'
            for keynum, key in enumerate(linedict.keys()):
                if keynum != 0:
                    convertline += ', '
                if key == 'data':
                    convertline += 'element data:'
                    for ele in linedict[key]:
                        convertline += ('\t' + _np.str(ele))
                else:
                    convertline += (key + ': ' + _np.str(linedict[key]))
                if keynum == len(linedict.keys()):
                    convertline += '.'
            self.Writer.DebugPrintout(convertline)

            if self.Transport.convprops.combineDrifts:
                if lastElementWasADrift and linedict['elementnum'] != 3.0 and linedict['elementnum'] < 20.0:
                    # write possibly combined drift
                    self.Writer.DebugPrintout('\n\tConvert delayed drift(s)')
                    self.Drift(linedictDrift)
                    lastElementWasADrift = False
                    self.Writer.DebugPrintout('\n\tNow convert element number' + _np.str(linenum))

            if linedict['elementnum'] == 15.0:
                self.UnitChange(linedict)
            if linedict['elementnum'] == 20.0:
                self.ChangeBend(linedict)
            if linedict['elementnum'] == 1.0:  # Add beam on first definition
                if not self.Transport.convprops.beamdefined:
                    self.DefineBeam(linedict)
                elif not self.Transport.convprops.dontSplit:  # Only update beyond first definition if splitting is permitted
                    self.DefineBeam(linedict)
            if linedict['elementnum'] == 3.0:
                if skipNextDrift:
                    skipNextDrift = False
                    continue
                if self.Transport.convprops.combineDrifts:
                    self.Writer.DebugPrintout('\tDelay drift')
                    if lastElementWasADrift:
                        linedictDrift['length'] += linedict['length']  # update linedictDrift
                        if not linedictDrift['name']:
                            linedictDrift['name'] = linedict['name']  # keep first non-empty name
                    else:
                        linedictDrift = linedict  # new linedictDrift
                        lastElementWasADrift = True
                else:
                    self.Drift(linedict)
            if linedict['elementnum'] == 4.0:
                self.Dipole(linedict)
            if linedict['elementnum'] == 5.0:
                self.Quadrupole(linedict)
            if linedict['elementnum'] == 6.0:
                if not self.Transport.convprops.typeCode6IsTransUpdate:
                    self.Collimator(linedict)
                    # Length gotten from next drift
                    if linedict['length'] > 0.0:
                        skipNextDrift = True
                else:
                    self.TransformUpdate(linedict)
            if linedict['elementnum'] == 12.0:
                self.Correction(linedict)
            if linedict['elementnum'] == 11.0:
                self.Acceleration(linedict)
            if linedict['elementnum'] == 13.0:
                self.Printline(linedict)
            if linedict['elementnum'] == 16.0:
                self.SpecialInput(linedict)
            if linedict['elementnum'] == 18.0:
                self.Sextupole(linedict)
            if linedict['elementnum'] == 19.0:
                self.Solenoid(linedict)

            # 9.  : 'Repetition' - for nesting elements
            if linedict['elementnum'] == 9.0:
                self.Writer.DebugPrintout('\tWARNING Repetition Element not implemented in converter!')
            if linedict['elementnum'] == 2.0:
                errorline = '\tLine is a poleface rotation which is handled by the previous or next dipole element.'
                self.Writer.DebugPrintout(errorline)

            self.Writer.DebugPrintout('\n')

        # OTHER TYPES WHICH CAN BE IGNORED:
        # 6.0.X : Update RX matrix used in TRANSPORT
        # 7.  : 'Shift beam centroid'
        # 8.  : Magnet alignment tolerances
        # 10. : Fitting constraint
        # 14. : Arbitrary transformation of TRANSPORT matrix
        # 22. : Space charge element
        # 23. : RF Cavity (Buncher), changes bunch energy spread

        # Write also last drift
        if self.Transport.convprops.combineDrifts:
            if lastElementWasADrift:
                self.Writer.DebugPrintout('\tConvert delayed drift(s)')
                self.Drift(linedictDrift)

    def _UpdateElementsFromFits(self):
        # Functions that update the elements in the element registry.
        # For debugging purposes, they return dictionaries of the element type,
        # length change details, and which parameters were updated and the values in a list which
        # follows the pattern of [parameter name (e.g. 'field'),oldvalue,newvalue]

        # Length update common to nearly all elements, seperate function to prevent duplication
        # fitIndex was used in the past. Pass in anyway in case of future need, just delete for now.
        def _updateLength(eleIndex, fitIndex, element):
            del fitIndex
            oldlength = self.Transport.ElementRegistry.elements[eleIndex]['length']
            lengthDiff = self.Transport.ElementRegistry.elements[eleIndex]['length'] - element['length']
            self.Transport.ElementRegistry.elements[eleIndex]['length'] = element['length']  # Update length
            self.Transport.ElementRegistry.length[eleIndex:] += lengthDiff  # Update running length of subsequent elements.
            self.Transport.ElementRegistry._totalLength += lengthDiff  # Update total length
            lendict = {'old': _np.round(oldlength, 5),
                       'new': _np.round(element['length'], 5)}
            return lendict

        def _updateDrift(eleIndex, fitIndex, element):
            eleDict = {'updated': False,
                       'element': 'Drift',
                       'params': []}

            # Only length can be varied
            if self.Transport.ElementRegistry.elements[eleIndex]['length'] != element['length']:
                lendict = _updateLength(eleIndex, fitIndex, element)
                eleDict['updated'] = True
                eleDict['length'] = lendict
            return eleDict

        def _updateQuad(eleIndex, fitIndex, element):
            eleDict = {'updated': False,
                       'element': 'Quadrupole',
                       'params': []}

            if self.Transport.ElementRegistry.elements[eleIndex]['data'][1] != element['data'][1]:  # Field
                oldvalue = self.Transport.ElementRegistry.elements[eleIndex]['data'][1]
                self.Transport.ElementRegistry.elements[eleIndex]['data'][1] = element['data'][1]
                eleDict['updated'] = True
                data = ['field', oldvalue, element['data'][1]]
                eleDict['params'].append(data)

            if self.Transport.ElementRegistry.elements[eleIndex]['length'] != element['length']:
                self.Transport.ElementRegistry.elements[eleIndex]['data'][0] = element['data'][0]  # Length in data
                lendict = _updateLength(eleIndex, fitIndex, element)
                eleDict['updated'] = True
                eleDict['length'] = lendict
            return eledict

        def _updateDipole(eleIndex, fitIndex, element):
            eleDict = {'updated': False,
                       'element': 'Dipole',
                       'params': []}

            # TODO: Need code in here to handle variation in poleface rotation. Not urgent for now.
            if self.Transport.ElementRegistry.elements[eleIndex]['data'][1] != element['data'][1]:  # Field
                oldvalue = self.Transport.ElementRegistry.elements[eleIndex]['data'][1]
                self.Transport.ElementRegistry.elements[eleIndex]['data'][1] = element['data'][1]
                eleDict['updated'] = True
                if self.Transport.machineprops.benddef:  # Transport can switch dipole input definition
                    par = 'field'
                else:
                    par = 'angle'
                data = [par, oldvalue, element['data'][3]]
                eleDict['params'].append(data)
            if self.Transport.ElementRegistry.elements[eleIndex]['length'] != element['length']:
                self.Transport.ElementRegistry.elements[eleIndex]['data'][0] = element['data'][0]  # Length in data
                lendict = _updateLength(eleIndex, fitIndex, element)
                eleDict['updated'] = True
                eleDict['length'] = lendict
            return eledict

        for index, name in enumerate(self.Transport.FitRegistry._uniquenames):
            fitstart = self.Transport.FitRegistry.GetElementStartSPosition(name)
            elestart = self.Transport.ElementRegistry.GetElementStartSPosition(name)
            fitindex = self.Transport.FitRegistry.GetElementIndex(name)
            eleindex = self.Transport.ElementRegistry.GetElementIndex(name)
            for fitnum, fit in enumerate(fitstart):
                for elenum, ele in enumerate(elestart):
                    if (_np.round(ele, 5) == _np.round(fit, 5)) and \
                            (not self.Transport.ElementRegistry.elements[eleindex[elenum]]['isZeroLength']):
                        fitelement = self.Transport.FitRegistry.elements[fitindex[fitnum]]
                        if fitelement['elementnum'] == 3:
                            eledict = _updateDrift(eleindex[elenum], fitindex[fitnum], fitelement)
                        elif fitelement['elementnum'] == 4:
                            eledict = _updateDipole(eleindex[elenum], fitindex[fitnum], fitelement)
                        elif fitelement['elementnum'] == 5:
                            eledict = _updateQuad(eleindex[elenum], fitindex[fitnum], fitelement)

                        if eledict['updated']:
                            self.Writer.DebugPrintout(
                                "\tElement " + _np.str(eleindex[elenum]) + " was updated from fitting.")
                            self.Writer.DebugPrintout("\tOptics Output line:")
                            self.Writer.DebugPrintout(
                                "\t\t'" + self.Transport.FitRegistry.lines[fitindex[fitnum]] + "'")
                            if eledict.has_key('length'):
                                lenline = "\t" + eledict['element'] + " length updated to "
                                lenline += _np.str(eledict['length']['new'])
                                lenline += " (from " + _np.str(eledict['length']['old']) + ")."
                                self.Writer.DebugPrintout(lenline)
                            for param in eledict['params']:
                                parline = "\t" + eledict['element'] + " " + param[0]
                                parline += " updated to " + _np.str(param[2]) + " (from " + _np.str(param[1]) + ")."
                                self.Writer.DebugPrintout(parline)
                                # self.Writer.DebugPrintout("\n")

                        break
