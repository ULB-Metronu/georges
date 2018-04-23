# pytransport.Reader - Transport input and output file readers.
# Version 1.0
# W. Shields and J. Snuverink
# william.shields.2010@live.rhul.ac.uk

"""
Reader

Readers for loading Transport input files and output files.
Can extract the individual lattice, optics, and fitting sections.

Classes:
Reader - a list of data read from Transport files.
ConversionData - a class for holding data during conversion.

"""

import numpy as _np
from Data import BDSData as _BDA
import _General


class Reader:
    """
    Class for reading Transport input and output files.
    The output file can be standard output or Beam output.

    """
    def __init__(self):
        self._allowedIndicatorLines = ['0  100', '0    0']
        self.optics = _Optics()

    def GetOptics(self, inputFile, inputType=None):
        """
        Extract the optics from a Transport output file.
        """
        if isinstance(inputType, _np.str):
            if inputType == 'beam':
                transdata = self.optics._getBeamOptics(inputFile)
                return transdata
            elif inputType == 'standard':
                transdata = self.optics._getStandardOptics(inputFile)
                return transdata
    
        f = open(inputFile)
        transdata = None
        for line in f:
            # remove any carriage returns (both Mac and Unix)
            line = line.rstrip('\r\n')
            splitline = _remove_blanks(line.split(' '))
            if splitline and splitline[0] == '*BEAM*':  # Is beam output
                # print "'*BEAM*' found in line " + _np.str(i+1)
                f.close()
                transdata = self.optics._getBeamOptics(inputFile)
                break
            elif line in self._allowedIndicatorLines:
                # print "'0    0' found in line " + _np.str(i+1)
                f.close()
                transdata = self.optics._getStandardOptics(inputFile)
                break
            else:
                pass
        if transdata is None:
            errorstring = "Could not find an indicator in the file for either a beam output file\n"
            errorstring += "(indicator = '*BEAM*), or a standard output file (indicator = '0    0').\n"
            errorstring += "Please check the input file or specify the input type with the type argument \n"
            errorstring += "in the get_output function. Note that the only accepted values for type are \n"
            errorstring += "'standard' or 'beam'."
            raise IOError(errorstring)
        return transdata

    def GetLattice(self, inputFile):
        """
        Function to extract the lattice from a standard output file.
        """
        flist = _LoadFile(inputFile)
        latticestart = 0
        latticeend = 0
        foundlatticestart = False
        foundlatticeend = False
        lattice = ['OUTPUT LATTICE']
        for linenum, line in enumerate(flist):
            if line in self._allowedIndicatorLines:
                if not foundlatticestart:
                    latticestart = linenum+1
                    if flist[linenum+1] == '0INDICATOR VALUE WRONG OR MISSING - ZERO ASSUMED':
                        latticestart += 1
                    foundlatticestart = True
            if line == '0SENTINEL':
                if not foundlatticeend:
                    latticeend = linenum
                    foundlatticeend = True
        if not foundlatticestart:
            if not foundlatticeend:
                raise IOError('No lattice found in ' + inputFile + '.')
            else:
                errorstring = 'The end of a lattice (line = "0SENTINEL") was found at line ' + _np.str(latticeend + 1) + ',\n'
                errorstring += 'but the start of a lattice (line = "0    0") was not found. Please check the input file.'
                raise IOError(errorstring)
        elif not foundlatticeend:
                errorstring = 'The start of a lattice (line = "0    0") was found at line ' + _np.str(latticestart - 1) + ',\n'
                errorstring += 'but the end of a lattice (line = "0SENTINEL") was not found. Please check the input file.'
                raise IOError(errorstring)
        else:
            lattice.extend(flist[latticestart:latticeend])
        return lattice

    def GetFits(self, inputFile):
        """
        Function to get the fit routine data from the standard transport output.
        Returns two lists, the first with the direct output from the fitting data,
        the second with the first line of each element in the output data, which contains the
        element parameters with their fitted values.
        """
        flist = _LoadFile(inputFile)
        fitstart = 0
        fitend = 0
        foundfitstart = False
        foundfitend = False
        fits = []
        for linenum, line in enumerate(flist):
            if not foundfitstart:
                if line == '0SENTINEL':
                    fitstart = linenum
                    foundfitstart = True
            if not foundfitend:
                splitline = _remove_blanks(line.split(' '))
                if splitline and splitline[0] == '*BEAM*' and not foundfitend:
                    fitend = linenum
                    foundfitend = True
        if not foundfitstart:
            print('No fitting output found.')
            return None
        elif foundfitstart and not foundfitend:
                errorstring = 'The start of the fitting output (first line containing "0SENTINEL") was found at line ' + _np.str(fitstart-1) + ',\n'
                errorstring += 'but the end of the fitting output (first line containing "*BEAM*") was not found. Please check the input file.'
                raise IOError(errorstring)
        fits.extend(flist[fitstart:fitend])

        output = self.optics._getOptics(flist, inputFile)
        fitres = [element[0] for element in output]

        return fits, fitres

    def GetLatticeAndOptics(self, inputFile):
        """
        Function to extract the lattice and optics from a standard output file.
        """
        flist = _LoadFile(inputFile)
        lattice = self.GetLattice(inputFile)
        optics = self.optics._getOptics(flist, inputFile)
        return lattice, optics
    

class _Optics:
    """
    Class for reading optics from Transport output files.
    The optics can be from standard output or Beam output.
    """
    def __init__(self):
        self.transdata = {
            'Sigma_x'   : [],
            'Sigma_xp'  : [],
            'Sigma_y'   : [],
            'Sigma_yp'  : [],
            'S'         : [],
            'Alpha_x'   : [],
            'Alpha_y'   : [],
            'Beta_x'    : [],
            'Beta_y'    : [],
            'Emitt_x'   : [],
            'Emitt_y'   : [],
            'Disp_x'    : [],
            'Disp_y'    : [],
            'Sigma_p'   : [],
            'Momentum'  : [],
            'E'         : [],  # kinetic energy
            'Name'      : [],
            'Type'      : []
            }
        # TODO: some unit for now, needs to be extracted from output / convert.py !
        self.transunits = {
            'Sigma_x'   : 'mm',
            'Sigma_xp'  : 'mrad',
            'Sigma_y'   : 'mm',
            'Sigma_yp'  : 'mrad',
            'S'         : 'm',
            'Alpha_x'   : '',
            'Alpha_y'   : '',
            'Beta_x'    : 'mm / mrad',
            'Beta_y'    : 'mm / mrad',
            'Emitt_x'   : 'mm mrad',
            'Emitt_y'   : 'mm mrad',
            'Disp_x'    : '',
            'Disp_y'    : '',
            'Sigma_p'   : 'MeV/c',
            'Momentum'  : 'MeV/c',
            'E'         : 'MeV',  # kinetic energy
            'Name'      : '',
            'Type'      : ''
            }

    def _processBeamOptics(self, flist):
        """
        Process the optics from a Beam output file.
        """
        transdata = {
            'Sigma_x'   : [],
            'Sigma_xp'  : [],
            'Sigma_y'   : [],
            'Sigma_yp'  : [],
            'S'         : [],
            'Alpha_x'   : [],
            'Alpha_y'   : [],
            'Beta_x'    : [],
            'Beta_y'    : [],
            'Emitt_x'   : [],
            'Emitt_y'   : [],
            'Disp_x'    : [],
            'Disp_y'    : [],
            'Sigma_p'   : [],
            'Name'      : [],
            }
        num_elements = 0
        for elenum, element in enumerate(flist):
            if element == '':  # The first line of the section should be a blank line.
                try:
                    section = flist[elenum+1:elenum+13]

                    # Get the element name, type and start position (in s)
                    line0 = section[0]
                    line0split = line0.split('*')
                    # typestr = line0split[1]
                    reststr = line0split[2]
                    restsplit = reststr.split(' ')
                    restsplit = _remove_blanks(restsplit)
                    transdata['S'].append(_np.float(restsplit[2]))
                    try :
                        namestr = restsplit[4]
                    except IndexError:
                        namestr = ''
                    transdata['Name'].append(namestr)
                    # transdata['type'].append(typestr)

                    # Get sigma_x and sigma_xp
                    line3 = section[3].split(' ')
                    line3 = _remove_blanks(line3)
                    transdata['Sigma_x'].append(_np.float(line3[3])/1000)
                    transdata['Sigma_xp'].append(_np.float(line3[5])/1000)

                    # Get sigma_y and sigma_yp
                    line4 = section[4].split(' ')
                    line4 = _remove_blanks(line4)
                    transdata['Sigma_y'].append(_np.float(line4[3])/1000)
                    transdata['Sigma_yp'].append(_np.float(line4[5])/1000)
                    
                    # Get momentum spread
                    line5 = section[5].split(' ')
                    transdata['Sigma_p'].append(_np.float(line5[5])/100)
                    
                    # Get alfa and beta twiss transdata for x and y.
                    line7 = section[7].split(' ')
                    line7 = _remove_blanks(line7)
                    transdata['Alpha_x'].append(_np.float(line7[0]))
                    transdata['Beta_x'].append(_np.float(line7[1]))
                    transdata['Alpha_y'].append(_np.float(line7[3]))
                    transdata['Beta_y'].append(_np.float(line7[4]))

                    # Get horizontal and vertical dispersion
                    line10 = section[10].split(' ')
                    line10 = _remove_blanks(line10)
                    transdata['Disp_x'].append(_np.float(line10[2])/10)
                    transdata['Disp_y'].append(_np.float(line10[5])/10)
                    
                    # Terms for calculating the emittance.
                    term1x = _np.float(line3[3])**2
                    term2x = (_np.float(line10[2])*(_np.float(line5[5])/100))**2
                    term1y = _np.float(line4[3])**2
                    term2y = (_np.float(line10[5])*(_np.float(line5[5])/100))**2

                    # Get horizontal and vertical emittance
                    emittx = (term1x - term2x) / _np.float(line7[1])
                    emitty = (term1y - term2y) / _np.float(line7[4])
                    transdata['Emitt_x'].append(emittx)
                    transdata['Emitt_y'].append(emitty)
                
                    num_elements += 1
                except ValueError:
                    errstr = "Could not process section beginning at line " + _np.str(elenum) + " : "
                    print(errstr)
                    print " "
                    print(element)
            elif element == "EOF -- rewind file":
                break

        data = _BDA()  # Now convert the dict into BDSData instance for final output.
        for keyName in transdata.keys():
            data._AddProperty(keyName)
        for i in range(num_elements):
            data.append(_GetElementData(i, transdata))

        return data

    def _processStandardOptics(self, elementlist, filename):
        """
        Process the optics from a standard output file.
        """
        if _General.CheckSingleLineOutputApplied(filename):
            optics = self._processStandardOpticsSingleLine(elementlist)
        else:
            optics = self._processStandardOpticsMultiLines(elementlist)
        return optics

    def _processStandardOpticsMultiLines(self, elementlist):
        """
        Process the optics from a standard output file when written to multiple lines.
        """
        # okElements=['BEAM','CORR','DRIFT','QUAD','SLIT','ADD TO BEAM','BEND','ROTAT','Z RO']
        notokElements = ['AXIS SHIFT']
        
        num_elements = 0
        # initialise momentum/energy since not given for every element
        momentum = 0.0
        energy = 0.0
        proton_mass = 938.272
        for element in elementlist:
            if (not isinstance(element, _np.str)) and (len(element) > 1):  # I.e not a fit or matrix-modifying element
                # type is in between * can have a space (for space charge *SP CH*)
                elementType = element[0].split('*')[1]
                if elementType not in notokElements:
                    # rest of first line split with spaces
                    splitline = _remove_blanks(element[0].split('*')[2].split(' '))
                    elename = splitline[1].strip('"')
                    if elementType == "BEAM" or elementType == "ACC":
                        momentum = _np.float(splitline[-2])
                        energy = _np.sqrt(proton_mass*proton_mass + momentum*momentum) - proton_mass
                    s       = _np.float(_remove_blanks(element[1].split(' '))[0])
                    sigx    = _np.float(_remove_blanks(element[1].split(' '))[3])
                    sigxp   = _np.float(_remove_blanks(element[2].split(' '))[1])
                    sigy    = _np.float(_remove_blanks(element[3].split(' '))[1])
                    sigyp   = _np.float(_remove_blanks(element[4].split(' '))[1])
                    sigt    = _np.float(_remove_blanks(element[5].split(' '))[1])
                    sigp    = _np.float(_remove_blanks(element[6].split(' '))[1])
                    r21     = _np.float(_remove_blanks(element[2].split(' '))[3])
                    r43     = _np.float(_remove_blanks(element[4].split(' '))[5])

                    dx = _GetTransformLineElements(element[8])[5]
                    dy = _GetTransformLineElements(element[10])[5]

                    self._SetTransportData(sigx, sigxp, sigy, sigyp, s, dx, dy, sigp, momentum, energy, elename,
                                           elementType, r21, r43)
                    num_elements += 1 

        data = _BDA()      # Now convert the dict into BDSData instance for final output.
        for keyName, unit in self.transunits.iteritems():
            data._AddProperty(keyName, unit)
        for i in range(num_elements):
            data.append(_GetElementData(i, self.transdata))

        return data

    def _SetTransportData(self, sigx, sigxp, sigy, sigyp, s, dx, dy, sigp, momentum, energy, elename, elementType,
                          r21, r43):
        """
        Set the beam data.
        """
        # Add/Subtract small amount if sin of phase space ellipse rotation is +/-one.
        # This comes from the output annoyingly rounding the code to one ,
        # which produces a div by zero later in the beta and gamma calculations.
        if r21 == 1.0:
            r21 -= 1e-4
        if r21 == -1.0:
            r21 += 1e-4

        if r43 == 1.0:
            r43 -= 1e-4
        if r43 == -1.0:
            r43 += 1e-4

        # Calculate twiss parameters
        xpint = _np.sqrt(sigxp ** 2 * (1 - r21 ** 2))
        ypint = _np.sqrt(sigyp ** 2 * (1 - r43 ** 2))

        ex = sigx * xpint
        ey = sigy * ypint

        if ex is 0:
            betx = 0
            gammax = 0
        else:
            betx = (sigx ** 2.0) / ex
            gammax = (sigxp ** 2.0) / ex

        if ey is 0:
            bety = 0
            gammay = 0
        else:
            bety = (sigy ** 2.0) / ey
            gammay = (sigyp ** 2.0) / ey

        if _np.isnan(betx):
            betx = 0
        if _np.isnan(bety):
            bety = 0
        if _np.isnan(gammax):
            gammax = 0
        if _np.isnan(gammay):
            gammay = 0

        alfx2 = (gammax * betx) - 1.0
        alfy2 = (gammay * bety) - 1.0
        if (alfx2 < 0):
            alfx2 = 0
        if (alfy2 < 0):
            alfy2 = 0

        alfx = _np.sqrt(alfx2)
        alfy = _np.sqrt(alfy2)

        self.transdata['Sigma_x'].append(sigx / 1000)  # convert to m
        self.transdata['Sigma_xp'].append(sigxp / 1000)  # convert to rad
        self.transdata['Sigma_y'].append(sigy / 1000)  # convert to m
        self.transdata['Sigma_yp'].append(sigyp / 1000)  # convert to rad
        self.transdata['S'].append(s)
        self.transdata['Alpha_x'].append(alfx)
        self.transdata['Alpha_y'].append(alfy)
        self.transdata['Beta_x'].append(betx)
        self.transdata['Beta_y'].append(bety)
        self.transdata['Emitt_x'].append(ex)
        self.transdata['Emitt_y'].append(ey)
        self.transdata['Disp_x'].append(dx)
        self.transdata['Disp_y'].append(dy)
        self.transdata['Sigma_p'].append(sigp)
        self.transdata['Momentum'].append(momentum)
        self.transdata['E'].append(energy)
        self.transdata['Name'].append(elename)
        self.transdata['Type'].append(elementType)

    def _processStandardOpticsSingleLine(self, elementlist):
        """
        Process the optics from a standard output file when written to single lines as specified
        by a 13. 19. element in Transport.
        """
        # seperate R matrix table from sigma matrix elements
        rMatrixElements = elementlist[-1]
        rMatrix = []

        # Second to last is column headers for R matrix table
        sMatrix = elementlist[:-2]

        for element in rMatrixElements[1:]:
            rMatrix.append(_remove_blanks(element.split(' ')))

        num_elements = 0
        momentum = 0.0
        energy = 0.0
        proton_mass = 938.272
        notokElements = ['AXIS SHIFT']
        okRElements = [3, 4, 5]  # ok element types for R matrix matching

        for element in sMatrix:
            if len(element) > 1:  # I.e not a fit or matrix-modifying element
                elementLine = _remove_blanks(element[0].split(' '))
                elementLine = _updateElementLine(elementLine)
                elementType = elementLine[0].strip('*')  # element type
                # typenum = _np.float(elementLine[1])

                # Get line with sigma data.
                # Element can have an additional line with fitting vary code data
                if len(element) == 3:
                    sigmaLine = element[2]
                else:
                    sigmaLine = element[1]

                if elementType not in notokElements:
                    elename = _removeIllegals(elementLine[2])  # remove illegal characters

                    if elementType == "BEAM" or elementType == "ACC":
                        momentum = _np.float(elementLine[-2])
                        energy = _np.sqrt(proton_mass*proton_mass + momentum*momentum) - proton_mass

                    if len(element) > 6:  # In case beam is defined before output format change.
                        s     = _np.float(_remove_blanks(element[1].split(' '))[0])
                        sigx  = _np.float(_remove_blanks(element[1].split(' '))[3])
                        sigxp = _np.float(_remove_blanks(element[2].split(' '))[1])
                        sigy  = _np.float(_remove_blanks(element[3].split(' '))[1])
                        sigyp = _np.float(_remove_blanks(element[4].split(' '))[1])
                        sigt  = _np.float(_remove_blanks(element[5].split(' '))[1])
                        sigp  = _np.float(_remove_blanks(element[6].split(' '))[1])

                        try:
                            r21 = _np.float(_remove_blanks(element[2].split(' '))[3])
                        except IndexError:
                            r21 = 0
                        try:
                            r43 = _np.float(_remove_blanks(element[4].split(' '))[5])
                        except IndexError:
                            r43 = 0
                    else:
                        s     = _np.float(_remove_blanks(sigmaLine.split(' '))[0])
                        sigx  = _np.float(_remove_blanks(sigmaLine.split(' '))[2])
                        sigxp = _np.float(_remove_blanks(sigmaLine.split(' '))[4])
                        sigy  = _np.float(_remove_blanks(sigmaLine.split(' '))[6])
                        sigyp = _np.float(_remove_blanks(sigmaLine.split(' '))[8])
                        sigt  = _np.float(_remove_blanks(sigmaLine.split(' '))[10])
                        sigp  = _np.float(_remove_blanks(sigmaLine.split(' '))[12])
                        r21   = _np.float(_remove_blanks(sigmaLine.split(' '))[14])
                        r43   = _np.float(_remove_blanks(sigmaLine.split(' '))[15])

                    dx = 0
                    dy = 0

                    # Find matching R matrix element and get dispersion
                    for rElement in rMatrix:
                        if _np.float(rElement[1]) in okRElements and (_np.float(rElement[0]) == s) \
                                and (rElement[2] == elename):
                            # Dispersion position dependent on existence of field strength in output
                            # Field strength written before first *, which should be the 4th element
                            if rElement.index('*') == 4:
                                dx = _np.float(rElement[15])
                                dy = _np.float(rElement[17])
                            else:
                                dx = _np.float(rElement[14])
                                dy = _np.float(rElement[16])

                    self._SetTransportData(sigx, sigxp, sigy, sigyp, s, dx, dy, sigp, momentum, energy, elename,
                                           elementType, r21, r43)
                    num_elements += 1

        data = _BDA()      # Now convert the dict into BDSData instance for final output.
        for keyName, unit in self.transunits.iteritems():
            data._AddProperty(keyName, unit)
        for i in range(num_elements):
            data.append(_GetElementData(i, self.transdata))
        
        return data

    def _getOptics(self, flist, filename):
        """
        Function to extract the output from a standard output file. The output will be an list of the lines
        for each element which contains the beam data. Each element should contain the R and TRANSPORT matrices
        which are necessary so the beam info can be calculated.
        """
        foundOpticsStart = False
        foundOpticsEnd = False
        foundRMatrixElementStart = False

        for linenum, line in enumerate(flist):
            splitline = _remove_blanks(line.split(' '))
            if splitline and splitline[0] == '*BEAM*' and not foundOpticsStart:
                opticsStart = linenum
                foundOpticsStart = True
            if splitline and splitline[0] == '0*LENGTH*':
                opticsEnd = linenum
                foundOpticsEnd = True
            if splitline and splitline[0] == '0POSITION':
                rMatrixElementStart = linenum
                foundRMatrixElementStart = True
    
        if not foundOpticsStart:
            if not foundOpticsEnd:
                raise IOError('No output found in ' + filename + '.')
            else:
                errorstring = 'The end of a lattice (line containing "0*LENGTH*") was found at ' \
                              'line ' + _np.str(opticsStart + 1)+',\n'
                errorstring += 'but the start of a lattice (first line containing "*BEAM*") was not found. ' \
                               'Please check the input file.'
                raise IOError(errorstring)
        elif not foundOpticsEnd:
                errorstring = 'The start of a lattice (first line containing "*BEAM*") was found at ' \
                              'line ' + _np.str(opticsStart - 1)+',\n'
                errorstring += 'but the end of a lattice (line containing "0*LENGTH*") was not found. ' \
                               'Please check the input file.'
                raise IOError(errorstring)
        else:
            output = flist[opticsStart:opticsEnd]

        # Append rest of the file which should only contain a table of R Matrix elements.
        if foundRMatrixElementStart:
            output.extend(flist[rMatrixElementStart:])

        # Split the list of all element data into their individual elements.
        elementlist = []
        elementstart = False
        elementend = False
        finalelement = False
        for linenum, line in enumerate(output):
            if linenum == (len(output)-1):
                finalelement = True
            try:
                if (line[1] == '*') or (line[:9] == '0POSITION'):
                    if line[2:11] == 'TRANSFORM':  # Is midway through element output
                        pass
                    elif elementstart is False:  # Current line must be start of the element
                        elementstart = True
                        elementstartline = linenum
                    elif elementstart is True:
                        elementend = True       # Otherwise the line must be the start of the next element
                        elementendline = linenum
                if elementstart and finalelement:
                    element = output[elementstartline:]
                    if element[-1][:2] == 'IO':
                        elementlist.append(element[:-1])
                        elementlist.append(element[-1])
                    else:
                        elementlist.append(element)
                if elementstart and elementend:  # If the start and end of the element are found, append and reset
                    element = output[elementstartline:elementendline]
                    if element[-1][:2] == 'IO':
                        elementlist.append(element[:-1])
                        elementlist.append(element[-1])
                    else:
                        elementlist.append(element)
                    elementend = False
                    elementstart = False
                if elementstart is False and elementend is False:  # Though if it's been reset, it must be because the
                    elementstart = True                            # current line is the start of next element, so set
                    elementstartline = linenum                     # the start line for the next element
            except IndexError:
                pass
        return elementlist

    def _getBeamOptics(self, inputFile):
        """
        Returns a BDSData instance of parameters from the input file.
        The input file is assumed to contain the beam data as output
        manually by TRANSPORT. As such, the structure of said output
        is assumed to remain unchanged from lattice to lattice. This code
        is designed to work with that output only. An example of that
        output would be:

        *DRIFT*      z = 34.372 m
        *SIGMA*
         Center:         0.000 mm    0.000 mrad    0.000 mm    0.000 mrad
         horz. Par. :   18.361 mm   23.134 mrad    0.997
         vert. Par. :    8.447 mm   10.242 mrad   -0.995
        *TWISS PARAMETERS* (for dp/p = 1.000 % )
           alfax:     betax:         alfay:     betay:
         -13.51800   10.75807 m      9.55361    7.92270 m
        *TRANSFORM 1*
            horz:                              vert:
          -14.22000   -0.64752   -9.32539      0.13646    0.97529    0.00000
          -18.40864   -0.90858  -10.33185     -1.18959   -1.17399    0.00000

        Anything other than this format will be read incorrectly or will
        not be read. There has to be a blank line between every element's
        output in the file.

        Note:
            It is assumed that the element (3,2) in the displayed *TRANSFORM 1*
            matrices corresponds to the dispersion (-9.32539 and 0.0000 in the above
            example). Some output however appears to have the correct magnitude, but
            incorrect sign. This doesn't affect the resulting beam size, but beware
            that a direct dispersion comparison to another lattice may appear incorrect.
        """
        flist = _LoadFile(inputFile)
        transdata = self._processBeamOptics(flist)
        return transdata

    def _getStandardOptics(self, inputFile):
        """
        Get the optics from a standard output file. Returns a pytransport.Data.BDSData object.
        """
        flist = _LoadFile(inputFile)
        optics = self._getOptics(flist, inputFile)
        transdata = self._processStandardOptics(optics, inputFile)
        return transdata


def _remove_blanks(line):
    """
    Removes any blanks from a string (blanks being '' and not white spaces).
    """
    newline = ''
    for element in line:
        if element != '':
            newline += element
            newline += ' '
    stripline = newline.split(' ')
    # remove last element as it will be blank (due to added space)
    return stripline[:-1]


def _removeIllegals(line):
    """
    Function to remove '' and stray characters from lines.
    """
    illegal = ['"', '', '(', ')']

    modLine = ''
    for element in line:
        if element not in illegal:
            modLine += element
    return modLine


def _split_negatives(line):
    newline = []
    for element in line:
        negpos = element.find('-')
        if negpos > 1:
            parts = element.split('-')
            newline.append(parts[0])
            newline.append('-' + parts[1])
        else:
            newline.append(element)
    return newline


def _LoadFile(inputfile):
    """
    Converts the input file into a list. The data has to be in a format other than
    handling line by line (using next()).

    The function for processing the file reads section by section rather than line
    by line. This may be an inefficient method but the input file should not be very
    large so it should not require a large amount of memory.
    """
    if inputfile == '':
        raise IOError('No file name supplied.')
    flist = []
    infile = open(inputfile)
    # Loop over lines and remove any carriage returns (both Mac and Unix)
    for line in infile:
        cleanline = line.rstrip('\r\n')
        flist.append(cleanline)
    infile.close()
    return flist


def _updateElementLine(line):
    if (line[0] == '*Z') and (line[1] == 'ROT*'):
        newline = []
        elementType = '*Z ROT*'
        newline.append(elementType)
        for element in line[2:]:
            newline.append(element)
    else:
        newline = line

    return newline


def _GetElementData(index, dataDict):
    # Function to get the data for each element, rather than each key.
    elementlist = [dataDict[keyName][index] for keyName in dataDict.keys()]
    return elementlist

def _GetTransformLineElements(line):
    elements = []
    for element in range(6):
        start = ((element + 1) * 10 + 1)
        end = start + 10
        eleVal = _remove_blanks(line[start:end].split(' '))[0]
        try:
            elements.append(_np.float(eleVal))
        except ValueError:
            elements.append(0)
    return elements


