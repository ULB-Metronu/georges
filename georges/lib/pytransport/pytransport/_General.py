# pytransport._General - general python scripts / tools
# Version 1.0
# W. Shields and J. Snuverink
# william.shields.2010@live.rhul.ac.uk

"""
General utilities for day to day housekeeping

Classes:
_Writer - a class used for writing any output during conversion.

"""
import numpy as _np
from scipy import constants as _con
import sys as _sys
import os as _os
import string as _string
import glob as _glob

import Reader as _Reader
from Data import _beamprops
from Data import ConversionData


class _Writer:
    """
    Class for writing terminal output, logfile output, and the converted machine. 
    Class is designed specifically for use with pytransport conversion. Class is 
    protected as it should be hidden.

    kwargs:
    debugOutput: bool, default = False. 
    If true, strings supplied to class functions for debug output will be written.
    writeToLog: bool, default = False. 
    If true, strings supplied to class functions will be written to a logfile.
    logfile: string, default = ''. Log file name.
    """
    def __init__(self, debugOutput=False, writeToLog=False, logfile=''):
        self.debug = debugOutput
        self.logfile = logfile
        self.outlog = writeToLog

    def Printout(self, line, outToTerminal=True):
        """
        Print line output string. Prints to output log if specified at class instantiation.
        Argument outToTerminal (bool, default = True) will print line to the terminal if True.
        """
        if outToTerminal:
            _sys.stdout.write(line+'\n')
        if self.outlog:
            if self.logfile == '':
                raise IOError("Invalid log file name: ''")
            logfile = open(self.logfile, 'a')
            logfile.write(line)
            logfile.write('\n')
            logfile.close()

    def DebugPrintout(self, line):
        """
        Print line debug output string to logfile only.
        """
        if self.debug:
            self.Printout(line, outToTerminal=False)

    def BeamDebugPrintout(self, beam, units):
        """
        Print the debug output strings for the beam definition.
        The beam must be a pytransport.Data._beamprops instance, and the units must be the
        pytransport.Data.ConversionData.units dict.
        """
        if not isinstance(beam, _beamprops):
            raise TypeError("Beam must be pytransport.Data._beamprops instance.")
        self.DebugPrintout('\tBeam definition :')
        self.DebugPrintout('\tdistrType = ' + beam.distrType)
        self.DebugPrintout('\tenergy = '  + _np.str(beam.tot_energy) + ' GeV')
        self.DebugPrintout('\tSigmaX = '  + _np.str(beam.SigmaX) + ' ' + units['x'])
        self.DebugPrintout('\tSigmaXP = ' + _np.str(beam.SigmaXP) + ' ' + units['xp'])
        self.DebugPrintout('\tSigmaY = '  + _np.str(beam.SigmaY) + ' ' + units['y'])
        self.DebugPrintout('\tSigmaYP = ' + _np.str(beam.SigmaYP) + ' ' + units['yp'])
        self.DebugPrintout('\tSigmaE = '  + _np.str(beam.SigmaE))
        self.DebugPrintout('\tSigmaT = '  + _np.str(beam.SigmaT))
        self.DebugPrintout('\t(Final brho = ' + _np.str(_np.round(beam.brho, 2)) + ' Tm)')
        self.DebugPrintout('\tTwiss Params:')
        self.DebugPrintout('\tBetaX = '  + _np.str(beam.betx) + ' ' + units['beta_func'])
        self.DebugPrintout('\tBetaY = '  + _np.str(beam.bety) + ' ' + units['beta_func'])
        self.DebugPrintout('\tAlphaX = ' + _np.str(beam.alfx))
        self.DebugPrintout('\tAlphaY = ' + _np.str(beam.alfy))
        self.DebugPrintout('\tEmittx = ' + _np.str(beam.emitx) + ' ' + units['emittance'])
        self.DebugPrintout('\tEmittY = ' + _np.str(beam.emity) + ' ' + units['emittance'])

    def ElementPrepDebugPrintout(self, elementType, numElements):
        """
        Print the debug output string as required in the element preparation stage.
        """
        debugString = "\tEntry is a " + elementType + ", adding to the element registry as element "
        debugString += numElements + "."
        self.DebugPrintout(debugString)

    def Write(self, convData, filename):
        """
        Write the converted TRANSPORT file to disk. A pytransport.Data.ConversionData 
        instance and filename (string)
        must be supplied.
        """
        if not isinstance(filename, _np.str):
            raise TypeError("Filename must be a string")
        if not isinstance(convData, ConversionData):
            raise TypeError("convData must be a pytransport.Data.ConversionData instance.")

        fname = ""
        directory = ""
        if convData.convprops.gmadoutput:
            fname = filename + '.gmad'
            directory = convData.convprops.gmadDir
        elif convData.convprops.madxoutput:
            fname = filename + '.madx'
            directory = convData.convprops.madxDir

        self.Printout('Writing to file: ' + directory + "/" + fname)
        if directory == "":
            convData.machine.Write(fname)
        else:
            if not CheckDirExists(directory):
                _os.mkdir(directory)
            _os.chdir(directory)
            convData.machine.Write(fname)
            _os.chdir('../')


def CheckDirExists(directory):
    dirs = _glob.glob('*/')
    if directory[-1] != '/':
        directory += '/'
    if directory in dirs:
        return True
    return False


def CheckIsAddition(line, filetype='input'):
    """
    Function to check if a BEAM line of TRANSPORT code is a beam definition or r.m.s addition.
    """
    # Output file is a standard format, any RMS addition line should always be 10 long.
    if filetype == 'output':
        if len(line) == 10:
            return True

    elif filetype == 'input':
        if len(line) > 8:
            if (line[8] == '0.') or (line[8] == '0'):
                return True
            else:
                return False
    else:
        raise ValueError("File type can only be input or output")


def CheckIsOutput(inputfile):
    """
    Function to check if a file is a standard TRANSPORT output file.
    Based upon existence of the lines::
        
         "0  XXX"
    
    being present, which represents the TRANSPORT indicator card line.
    X can be 0, 1, 2. Default is 0.
    """
    temp = _Reader.Reader()
    isOutput = False
    try:
        f = open(inputfile)
        for inputline in f:
            inputline = inputline.replace("\r", '')
            inputline = inputline.replace("\n", '')
            if inputline in temp._allowedIndicatorLines:
                isOutput = True
                break
        f.close()
    except IOError:
        raise IOError('Cannot open file.')
    return isOutput


def CheckIsSentinel(line):
    for element in line:
        if element[:8] == 'SENTINEL':
            return True
    return False


def CheckSingleLineOutputApplied(inputfile):
    """
    Function to check if the control element that print element output in
    a single line was successfully applied. Check needed as not all versions
    of TRANSPORT can run this type code.
    """
    reader = _Reader.Reader()
    flist = _Reader._LoadFile(inputfile)
    optics = reader.optics._getOptics(flist, inputfile)
    for element in optics:
        if element == 'IO: UNDEFINED TYPE CODE 13. 19. ;':
            return True
    return False


def ConvertBunchLength(transport, bunch_length):
    """
    Function to convert bunch length unit in TRANSPORT into seconds.
    """
    scale = transport.scale[transport.units['bunch_length'][0]]
    blmeters = bunch_length * scale  # Bunch length scaled to metres
    blseconds = blmeters / (transport.beamprops.beta*_con.c)  # Length converted to seconds
    return blseconds


def FindEndOfLine(line):  # Find the end of the line of code
    endpos = -1
    breakloop = False
    if isinstance(line, _np.str):
        for charnum, char in enumerate(line):
            if char == ';':
                endpos = charnum
                break
    elif isinstance(line, _np.ndarray):
        for index, ele in enumerate(line):
            for char in ele:
                if char == ';':
                    endpos = index
                    breakloop = True
                    break
            if breakloop:
                break
    return endpos


def GetComment(line):
    """
    Function to extract a comment from a line.
    Will only find one comment per line.
    """
    concat = ''
    for ele in line:
        concat += ele
        concat += ' '
    commstart = _string.find(concat, '(')
    commend = _string.find(concat, ')')
    if commstart != -1 and commend != -1:
        comment = concat[commstart+1: commend]
        gmadcomment = '! '+comment
    else:
        gmadcomment = None
    return gmadcomment


def GetElementData(line):
    data = []
    for index, ele in enumerate(line[1:]):
        if ele != '':
            try:
                data.append(_np.float(ele))
            except ValueError:
                pass
    return data


def GetFaceRotationAngles(data, linenum):

    def searchForAngle(linelist):
        angle = 0
        for line in linelist:
            try:
                elecode = _np.float(line[0])
                del elecode
            except ValueError:
                angle = 0
                break

            if _np.float(line[0]) == 4.0:
                break
            elif _np.float(line[0]) == 2.0:
                endof = FindEndOfLine(line[1])
                if endof != -1:
                    try:
                        angle = _np.round(_np.float(line[1][:endof]), 4)
                    except ValueError:
                        try:
                            angle = _np.round(_np.float(line[2][:endof]), 4)
                        except ValueError:
                            pass
                    else:
                        pass
                else:
                    try:
                        angle = _np.round(_np.float(line[1]), 4)
                    except ValueError:
                        try:
                            angle = _np.round(_np.float(line[2]), 4)
                        except ValueError:
                            pass
                    else:
                        pass
                break
            else:
                pass
        return angle

    lineList = [i for i in data[(linenum - 5):linenum]]
    lineList.reverse()  # Search for input poleface in reverse line order

    anglein = searchForAngle(lineList)
    angleout = searchForAngle(data[linenum + 1:(linenum + 6)])

    return anglein, angleout


def GetIndicator(data):
    """
    Function to read the indicator number. Must be 0, 1, or 2, where:
        0 is a new lattice
        1 is for fitting with the first lattice (from a 0 indicator file)
        2 if for a second fitting which suppresses the first fitting.
    """
    indc = 0
    linenum = 0
    for linenum, line in enumerate(data):
        if line[0] == '0':
            if line[1] == '\r' or '\n':
                indc = 0
                break
        if line[0] == '1':
            if line[1] == '\r' or '\n':
                indc = 1
                break
        if line[0] == '2':
            if line[1] == '\r' or '\n':
                indc = 2
                break
    return indc, linenum


def GetLabel(line):
    """
    Function to get element label from code line.
    """
    label = None
    for elenum, ele in enumerate(line):
        startslash = _string.find(ele, "/")
        startquote = _string.find(ele, "'")
        startequal = _string.find(ele, "=")
        startdbquote = _string.find(ele, '"')
        if startslash != -1:
            end = 1 + startslash + _string.find(ele[(startslash + 1):], "/")
            if end <= startslash:
                label = ele[startslash + 1:]
            else:
                label = ele[startslash + 1:end]
            break
        elif startquote != -1:
            end = 1 + startquote + _string.find(ele[(startslash + 1):], "'")
            if end <= startquote:
                label = ele[startquote + 1:]
            else:
                label = ele[startquote + 1:end]
            break
        elif startequal != -1:
            end = 1 + startequal + _string.find(ele[(startslash + 1):], "=")
            if end <= startequal:
                label = ele[startequal + 1:]
            else:
                label = ele[startequal + 1:end]
            break
        elif startdbquote != -1:
            end = 1 + startdbquote + _string.find(ele[(startdbquote + 1):], '"')
            if end <= startdbquote:
                label = ele[startdbquote + 1:]
            else:
                label = ele[startdbquote + 1:end]
            break
        else:
            label = None
    return label


def GetPreamble(data):  # Redundant until pybdsim can handle comments.
    """
    Function to read any preamble at the start of the TRANSPORT file.
    """
    # indc, linenum = GetIndicator(data)
    gmadpreamble = []
    # for line in self.Transport.data[:linenum-1]:
    for line in data:
        if line == '\r\n':
            pass
        else:
            gmadline = '!' + line
            gmadpreamble.append(gmadline)
    return gmadpreamble


def GetTypeNum(line):
    """
    Function to extract the element type number (type code).
    Written because element types can contain alphabetical
    characters when fits are used, e.g: 5.0A. Only the number
    is required, the use of fitting does not need to be known.
    """
    eleNum = line[0]
    characNum = 0
    if len(eleNum) > 2:
        for characNum in range(len(eleNum[2:])):
            try:
                converted = _np.float(eleNum[:characNum+2])
                del converted
            except ValueError:
                break
        typeNum = _np.float(eleNum[:characNum+2])
    else:
        typeNum = _np.float(eleNum)
    return typeNum


def JoinSplitLines(linenum, lattice):
    firstline = lattice[linenum].replace(';', '')
    latticeline = firstline  # Copy for later
    firstline = _np.array(firstline.split(' '), dtype=_np.str)
    firstline = RemoveIllegals(firstline)
    numericals = []
    nonnumericals = []
    # Keep entries that are strings of numbers
    for i in firstline:
        try:
            number = _np.float(i)
            numericals.append(_np.str(number))
        except ValueError:
            nonnumericals.append(i)

    # Number of numerical elements minus the first which should be the entry type number.
    # This is bascially a way of extracting any label or comments.
    numelements = len(numericals) - 1

    secline = lattice[linenum + 1].replace(';', '')
    secline = _np.array(secline.split(' '), dtype=_np.str)
    secline = RemoveIllegals(secline)
    secnumericals = []

    for i in secline:
        try:
            number = _np.float(i)
            secnumericals.append("%.4f" % number)
        except ValueError:
            pass

    # Second line should be 15 minus number of numerical elements from prev line.
    # This is done to skip erroneous numbers in the line such as '000' which have
    # appeared when lines have been split.
    secline = secnumericals[-15 + numelements:]
    numericals.extend(secline)

    # Add to latticeline so as to appear like one single line in the file
    seclinetxt = ""
    for i in secline:
        newline = "     " + i
        seclinetxt += newline
    latticeline += seclinetxt

    # Add name to output
    if len(nonnumericals) == 1:
        numericals.append(nonnumericals[0])
    line = _np.array(numericals)
    return latticeline, line


def OutputFitsToRegistry(transport, outputdata):
    isLegal = {'*DRIFT*': 3.0,
               '*QUAD*': 5.0,
               '*BEND*': 4.0}

    for line in outputdata:
        append = False
        linedict = {'elementnum': 0.0,
                    'name': '',
                    'length': 0.0}
        data = RemoveIllegals(line.split(' '))
        eledata = GetElementData(data)
        label = GetLabel(data)
        if data[0] in isLegal:
            linedict['elementnum'] = isLegal[data[0]]
            linedict['name'] = label
            linedict['data'] = eledata[1:]  # first value is elementnum.
            linedict['length'] = eledata[1]
            append = True

        # Only add an element with a name to the fitting registry.
        # (Element has to be named to be varied in the fitting routine).
        # Otherwise update the total length of the machine.
        if append and (label is not None) and (label != ''):
            transport.FitRegistry.AddToRegistry(linedict, line)
        else:
            transport.FitRegistry.UpdateLength(linedict)
    return transport


def ProcessFits(fits):  # redundant
    # First split the fitting output into its respective sections (input problem step).
    fitsections = []
    fitsstarts = []
    # Start line of each section
    for linenum, line in enumerate(fits):
        if line.startswith('1'):
            fitsstarts.append(linenum)

    for secnum in range(len(fitsstarts)):
        if secnum + 1 < len(fitsstarts):
            section = fits[fitsstarts[secnum]:fitsstarts[secnum + 1]]
        else:
            section = fits[fitsstarts[secnum]:]
        lines = []
        for line in section:
            lines.append(RemoveIllegals(line.split(' ')))
        fitsections.append(lines)

    magnetlines = []
    for section in fitsections:
        for line in section:
            if (len(line) > 0) and (line[0][0] == '*' and line[0][-1] == '*') and line[0] != '*FIT*':
                magnetlines.append(line)


def RemoveFileExt(inputfile):
    """
    Remove the file extension from the input file name. Only works on known extensions.
    """
    exts = [".txt", ".dat", ".DAT", ".TXT"]
    if inputfile[-4:] in exts:
        return inputfile[:-4]
    return inputfile


def RemoveIllegals(line):
    """
    Function to remove '' and stray characters from lines.
    """
    illegal = ['"', '', '(', ')']

    linelist = [element for element in line if element not in illegal]
    line = _np.array(linelist)
    return line


def RemoveLabel(line):
    """
    Function to remove the label from a line.
    """
    label, elenum = GetLabel(line)
    if label is not None:
        element = line[elenum]
        lablen = len(label)
        newval = element
        for index in range(len(element)):
            if element[index:index + lablen] == label:
                prelabel = element[:index - 1]
                postlabel = element[index + lablen + 1:]
                newval = prelabel + ' ' + postlabel
                break
        line[elenum] = newval
    return line


def RemoveSpaces(line):
    elementlist = []
    for i in line:
        if (i != '') and (i != ' '):
            elementlist.append(i)
    line = _np.array(elementlist)
    return line


def ScaleToMeters(transport, quantity):
    """
    Function to scale quantity (string) to meters, returns conversion factor.
    """
    if transport.units[quantity] != 'm':
        conversionFactor = transport.scale[transport.units[quantity][0]]
    else:
        conversionFactor = 1
    return conversionFactor


def UpdateEnergyFromMomentum(transport, momentum):
    """
    Function to calculate (from momentum):
        Total Energy
        Kinetic Energy
        Momentum
        Lorentz factor (gamma)
        Velocity (beta)
        Magnetic rigidity (brho)
    """
    momentum = _np.float(momentum)
    transport.beamprops.momentum = momentum
    p_mass = transport.beamprops.mass  # Particle rest mass (in GeV)
    scaling = 1
    mom_in_ev = momentum

    mom_unit = transport.units['p_egain']
    if mom_unit != 'eV':
        scaling = 1e9 / transport.scale[mom_unit[0]]     # Scaling relative to mom. unit
        mom_in_ev = momentum * transport.scale[mom_unit[0]]
    elif mom_unit == 'eV':
        scaling = 1e9                               # Scaling relative to 1 eV
        mom_in_ev = momentum
    p_mass *= scaling                               # Scale particle rest mass
    energy = _np.sqrt((p_mass**2) + (momentum**2))
    transport.beamprops.tot_energy = energy
    transport.beamprops.tot_energy_current = energy
    transport.beamprops.k_energy = energy - p_mass
    transport.beamprops.gamma = energy / p_mass
    transport.beamprops.beta = _np.sqrt((1.0 - (1.0 / transport.beamprops.gamma**2)))
    transport.beamprops.brho = mom_in_ev / _con.c
    return transport


def UpdateMomentumFromEnergy(transport, k_energy):
    """
    Function to calculate (from kinetic energy):
        Total Energy
        Kinetic Energy
        Momentum
        Lorentz factor (gamma)
        Velocity (beta)
        Magnetic rigidity (brho)
    """
    scaling = 1  # defaults
    mom_in_ev = 0
    k_energy = _np.float(k_energy)

    transport.beamprops.k_energy = k_energy
    p_mass = transport.beamprops.mass  # Particle rest mass (in GeV)

    e_unit = transport.units['p_egain']
    if e_unit != 'eV':
        scaling = 1e9 / transport.scale[e_unit[0]]  # Scaling relative to mom. unit
    elif e_unit == 'eV':
        scaling = 1e9  # Scaling relative to 1 eV
    p_mass *= scaling  # Scale particle rest mass

    # energy = _np.sqrt((p_mass**2 * _con.c**2) + (momentum**2 * _con.c**2)) / _con.c

    transport.beamprops.tot_energy_current = k_energy + p_mass
    transport.beamprops.momentum = _np.sqrt((transport.beamprops.tot_energy_current ** 2) - (p_mass ** 2))

    transport.beamprops.gamma = transport.beamprops.tot_energy_current / p_mass
    transport.beamprops.beta = _np.sqrt((1.0 - (1.0 / transport.beamprops.gamma ** 2)))

    if e_unit != 'eV':
        mom_in_ev = transport.beamprops.momentum * transport.scale[e_unit[0]]
    elif e_unit == 'eV':
        mom_in_ev = transport.beamprops.momentum

    transport.beamprops.brho = mom_in_ev / _con.c
    return transport
