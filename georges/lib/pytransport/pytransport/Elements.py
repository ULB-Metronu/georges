# pytransport.Elements - tools for Transport element conversion
# Version 1.0
# W. Shields and J. Snuverink
# william.shields.2010@live.rhul.ac.uk

"""
Functions to help convert elements from Transport to gmad/madx.
Code inherited in classes in pybdsim.Convert.Transport2Gmad and pymadx.Convert.Transport2Madx.

Classes:
Elements - a class used to convert different element types.

"""

import numpy as _np

import _General
from _General import _Writer
from Data import ConversionData as _convData


class Elements:
    """
    Class for converting different Transport element types into gmad/madx format.

    Required: transportData, dtype = pytransport.Data.ConversionData

    """
    def __init__(self, transportData):
        if not isinstance(transportData, _convData):
            raise TypeError("transportData must be a pytransport.Data.ConversionData instance")
        self.Transport = transportData
        logfileName = _General.RemoveFileExt(self.Transport.convprops.file) + '_conversion.log'
        self.Writer = _Writer(debugOutput=self.Transport.convprops.debug,
                              writeToLog=self.Transport.convprops.outlog,
                              logfile=logfileName)

    def DefineBeam(self, linedict):
        if linedict['isAddition']:
            if self.Transport.convprops.debug:
                self.Writer.DebugPrintout('\tIgnoring beam rms addition.')
            return
        if self.Transport.convprops.beamdefined and not self.Transport.convprops.dontSplit:
            self.Writer.Printout('Beam redefinition found. Writing previous section to file.')
            self.Writer.Printout('Splitting into multiple machines.')
            self.Transport.convprops.numberparts += 1
            self.Transport.AddBeam()
            if self.Transport.convprops.gmadoutput:
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
            self.Transport._NewMachines()
            self.Transport.convprops.correctedbeamdef = False

        momentum = linedict['momentum']

        self.Transport.convprops.beamdefined = True

        # Convert momentum to energy and set distribution params.
        self.Transport = _General.UpdateEnergyFromMomentum(self.Transport, momentum)
        self.Transport.beamprops.SigmaX  = _np.float(linedict['Sigmax'])
        self.Transport.beamprops.SigmaY  = _np.float(linedict['Sigmay'])
        self.Transport.beamprops.SigmaXP = _np.float(linedict['Sigmaxp'])
        self.Transport.beamprops.SigmaYP = _np.float(linedict['Sigmayp'])
        self.Transport.beamprops.SigmaE  = _np.float(linedict['SigmaE']) * 0.01 * (self.Transport.beamprops.beta**2)  # Convert from percentage mom spread to absolute espread
        self.Transport.beamprops.SigmaT  = _General.ConvertBunchLength(self.Transport, _np.float(linedict['SigmaT']))  # Get bunch length in seconds.

        # Calculate Initial Twiss params
        try:
            self.Transport.beamprops.betx = self.Transport.beamprops.SigmaX / self.Transport.beamprops.SigmaXP
        except ZeroDivisionError:
            self.Transport.beamprops.betx = 0
        try:
            self.Transport.beamprops.bety = self.Transport.beamprops.SigmaY / self.Transport.beamprops.SigmaYP
        except ZeroDivisionError:
            self.Transport.beamprops.bety = 0
        self.Transport.beamprops.emitx = self.Transport.beamprops.SigmaX * self.Transport.beamprops.SigmaXP / 1000.0
        self.Transport.beamprops.emity = self.Transport.beamprops.SigmaY * self.Transport.beamprops.SigmaYP / 1000.0

        self.Writer.BeamDebugPrintout(self.Transport.beamprops, self.Transport.units)

    def Drift(self, linedict):
        driftlen = linedict['length']
        if driftlen <= 0:
            self.Writer.DebugPrintout('\tZero or negative length element, ignoring.')
            return

        lenInM = driftlen * _General.ScaleToMeters(self.Transport, 'element_length')  # length in metres

        self.Transport.machineprops.drifts += 1
        elementid = ''
        if self.Transport.convprops.keepName:
            elementid = linedict['name']
        if not elementid:  # check on empty string
            elementid = 'DR'+_np.str(self.Transport.machineprops.drifts)

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddDrift(name=elementid, length=lenInM)

        self.Writer.DebugPrintout('\tConverted to:')
        self.Writer.DebugPrintout('\t' + 'Drift ' + elementid + ', length ' + _np.str(lenInM) + ' m')

    def Dipole(self, linedict):
        linenum = linedict['linenum']
        dipoledata = linedict['data']
        length = dipoledata[0]  # First two non-blanks must be the entries in a specific order.

        # Get poleface rotation
        e1 = linedict['e1'] * ((_np.pi / 180.0)*self.Transport.machineprops.bending)  # Entrance pole face rotation.
        e2 = linedict['e2'] * ((_np.pi / 180.0)*self.Transport.machineprops.bending)  # Exit pole face rotation.

        if e1 != 0:
            self.Writer.DebugPrintout('\tPreceding element (' + _np.str(linenum-1) + ') provides an entrance poleface rotation of ' + _np.str(_np.round(e1, 4)) + ' rad.')
        if e2 != 0:
            self.Writer.DebugPrintout('\tFollowing element (' + _np.str(linenum+1) + ') provides an exit poleface rotation of ' + _np.str(_np.round(e2, 4)) + ' rad.')

        # Fringe Field Integrals
        fintval = 0
        fintxval = 0
        if e1 != 0:
            fintval = self.Transport.machineprops.fringeIntegral
        if e2 != 0:
            fintxval = self.Transport.machineprops.fringeIntegral
        if (fintval != 0) or (fintxval != 0):
            self.Writer.DebugPrintout('\tA previous entry set the fringe field integral K1=' + _np.str(self.Transport.machineprops.fringeIntegral) + '.')
        if (fintval != 0) and (e1 != 0):
            self.Writer.DebugPrintout('\tSetting fint=' + _np.str(fintval) + '.')
        if (fintxval != 0) and (e2 != 0):
            self.Writer.DebugPrintout('\tSetting fintx=' + _np.str(fintxval) + '.')

        hgap = self.Transport.machineprops.dipoleVertAper * self.Transport.scale[self.Transport.units['bend_vert_gap'][0]] * 10

        # Calculate bending angle
        if self.Transport.machineprops.benddef:
            bfield = dipoledata[1]
            field_in_Gauss = bfield * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Scale to Gauss
            field_in_Tesla = field_in_Gauss * 1e-4                                  # Convert to Tesla
            if field_in_Tesla == 0:
                angle = 0                                                           # zero field = zero angle
            else:
                rho = self.Transport.beamprops.brho / (_np.float(field_in_Tesla))         # Calculate bending radius.
                angle = (_np.float(length) / rho) * self.Transport.machineprops.bending   # for direction of bend
            self.Writer.DebugPrintout('\tbfield = ' + _np.str(field_in_Gauss) + ' kG')
            self.Writer.DebugPrintout('\tbfield = ' + _np.str(field_in_Tesla) + ' T')
            self.Writer.DebugPrintout('\tCorresponds to angle of ' + _np.str(_np.round(angle, 4)) + ' rad.')
        else:
            angle_in_deg = dipoledata[1]
            angle = angle_in_deg * (_np.pi/180.) * self.Transport.machineprops.bending

        # Convert element length
        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        self.Transport.machineprops.dipoles += 1
        elementid = ''
        if self.Transport.convprops.keepName:
            elementid = linedict['name']
        if not elementid:  # check on empty string
            elementid = 'BM' + _np.str(self.Transport.machineprops.dipoles)

        # pybdsim and pymadx are the same. Check for non zero pole face rotation.
        if (e1 != 0) and (e2 != 0):
            self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM, angle=_np.round(angle, 4), e1=_np.round(e1, 4), e2=_np.round(e2, 4), fint=fintval, fintx=fintxval, hgap=hgap)
        elif (e1 != 0) and (e2 == 0):
            self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM, angle=_np.round(angle, 4), e1=_np.round(e1, 4), fint=fintval, fintx=0, hgap=hgap)
        elif (e1 == 0) and (e2 != 0):
            self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM, angle=_np.round(angle, 4), e2=_np.round(e2, 4), fint=0, fintx=fintxval, hgap=hgap)
        else:
            self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM, angle=_np.round(angle, 4))

        # Debug output
        if (e1 != 0) and (e2 != 0):
            polefacestr = ', e1= ' + _np.str(_np.round(e1, 4)) + ' rad, e2= ' + _np.str(_np.round(e2, 4)) + ' rad'
        elif (e1 != 0) and (e2 == 0):
            polefacestr = ', e1= ' + _np.str(_np.round(e1, 4)) + ' rad'
        elif (e1 == 0) and (e2 != 0):
            polefacestr = ', e2= ' + _np.str(_np.round(e2, 4)) + ' rad'
        else:
            polefacestr = ''

        if (fintval != 0) and (fintxval != 0):
            fringestr = ' , fint= ' + _np.str(fintval) + ', fintx= ' + _np.str(fintxval)
        elif (fintval != 0) and (fintxval == 0):
            fringestr = ' , fint= ' + _np.str(fintval)
        elif (fintval == 0) and (fintxval != 0):
            fringestr = ' , fintx= ' + _np.str(fintxval)
        else:
            fringestr = ''

        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Dipole ' + elementid + ', length= ' + _np.str(lenInM) + ' m, angle= ' + \
                      _np.str(_np.round(angle, 4)) + ' rad' + polefacestr + fringestr
        self.Writer.DebugPrintout('\t' + debugstring)

    def ChangeBend(self, linedict):
        """
        Function to change the direction of the dipole bend. Can be a direction other than horizontal (i.e != n*pi).
        """
        # NOT FULLY TESTED.
        angle = linedict['angle']
        rotation = False
        elementid = ''

        self.Transport.machineprops.angle = _np.float(angle)
        if self.Transport.machineprops.angle >= 360:
            self.Transport.machineprops.angle = _np.mod(self.Transport.machineprops.angle, 360)
        if self.Transport.machineprops.angle <= -360:
            self.Transport.machineprops.angle = _np.mod(self.Transport.machineprops.angle, -360)

        if self.Transport.machineprops.angle == 180 or self.Transport.machineprops.angle == -180:  # If 180 degrees, switch bending angle
            self.Transport.machineprops.bending *= -1

        elif self.Transport.machineprops.angle != 0:  # If not 180 degrees, use transform3d.
            # self.machineprops.angle *= -1
            # For conversion to correct direction. Eg in TRANSPORT -90 is upwards, in BDSIM, 90 is upwards.
            anginrad = self.Transport.machineprops.angle * (_np.pi / 180)
            self.Transport.machineprops.transforms += 1
            if self.Transport.convprops.keepName:
                elementid = linedict['name']
            if not elementid:  # check on empty string
                elementid = 't' + _np.str(self.Transport.machineprops.transforms)

            # only call for gmad, warning for madx
            if self.Transport.convprops.gmadoutput:
                self.Transport.machine.AddTransform3D(name=elementid, psi=anginrad)
            elif self.Transport.convprops.madxoutput:
                self.Writer.DebugPrintout('\tWarning, MadX Builder does not have Transform 3D!')

            rotation = True

        if rotation:
            self.Writer.DebugPrintout('\tConverted to:')
            debugstring = '\tTransform3D ' + elementid + ', angle ' + _np.str(_np.round(self.Transport.machineprops.angle, 4)) + ' rad'
            self.Writer.DebugPrintout('\t'+debugstring)
        elif self.Transport.machineprops.angle == 180:
            self.Writer.DebugPrintout('\tBending direction set to Right')
        elif self.Transport.machineprops.angle == -180:
            self.Writer.DebugPrintout('\tBending direction set to Left')

    def Quadrupole(self, linedict):
        quaddata = linedict['data']
        length = quaddata[0]        # First three non-blanks must be the entries in a specific order.
        field_at_tip = quaddata[1]  # Field in TRANSPORT units
        pipe_rad = quaddata[2]      # Pipe Radius In TRANSPORT units

        field_in_Gauss = field_at_tip * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  # Convert to Tesla

        pipe_in_metres = pipe_rad * _General.ScaleToMeters(self.Transport, 'bend_vert_gap')
        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        field_gradient = (field_in_Tesla / pipe_in_metres) / self.Transport.beamprops.brho  # K1 in correct units

        self.Transport.machineprops.quads += 1

        elementid = ''
        if self.Transport.convprops.keepName:
            elementid = linedict['name']
        if not elementid:  # check on empty string
            if field_gradient > 0:
                elementid = 'QF' + _np.str(self.Transport.machineprops.quads)
            elif field_gradient < 0:
                elementid = 'QD' + _np.str(self.Transport.machineprops.quads)
            else:
                elementid = 'NULLQUAD' + _np.str(self.Transport.machineprops.quads)  # For K1 = 0.

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddQuadrupole(name=elementid, length=lenInM, k1=_np.round(field_gradient, 4))

        string1 = '\tQuadrupole, field in gauss = ' + _np.str(field_in_Gauss) + ' G, field in Tesla = ' + _np.str(field_in_Tesla) + ' T.'
        string2 = '\tBeampipe radius = ' + _np.str(pipe_in_metres) + ' m. Field gradient = '+ _np.str(field_in_Tesla/pipe_in_metres) + ' T/m.'
        string3 = '\tBrho = ' + _np.str(_np.round(self.Transport.beamprops.brho, 4)) + ' Tm. K1 = ' +_np.str(_np.round(field_gradient, 4)) + ' m^-2'
        self.Writer.DebugPrintout(string1)
        self.Writer.DebugPrintout(string2)
        self.Writer.DebugPrintout(string3)
        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Quadrupole ' + elementid + ', length= ' + _np.str(lenInM) + ' m, k1= ' + _np.str(_np.round(field_gradient, 4)) + ' T/m'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Collimator(self, linedict):
        """
        A Function that writes the properties of a collimator element
        Only added for gmad, not for madx!
        """
        if linedict['length'] <= 0:
            self.Writer.DebugPrintout('\tZero or negative length element, ignoring.')
            return
        colldata = linedict['data']

        # Determine which entry is for horiz. and vert.
        aperx = self.Transport.machineprops.beampiperadius
        apery = self.Transport.machineprops.beampiperadius
        if _np.float(colldata[0]) == 1.0:
            aperx = colldata[1]
        elif _np.float(colldata[0]) == 3.0:
            apery = colldata[1]

        if len(colldata) > 2:
            if _np.float(colldata[2]) == 1.0:
                aperx = colldata[3]
            elif _np.float(colldata[2]) == 3.0:
                apery = colldata[3]
        aperx = _np.float(aperx)
        apery = _np.float(apery)

        lenInM = linedict['length'] * _General.ScaleToMeters(self.Transport, 'element_length')
        aperx_in_metres = aperx * _General.ScaleToMeters(self.Transport, 'x')
        apery_in_metres = apery * _General.ScaleToMeters(self.Transport, 'y')

        self.Transport.machineprops.collimators += 1
        elementid = ''
        if self.Transport.convprops.keepName:
            elementid = linedict['name']
        if not elementid:  # check on empty string
            elementid = 'COL'+_np.str(self.Transport.machineprops.collimators)

        collimatorMaterial = 'copper'  # Default in BDSIM, added to prevent warnings
        # only call for gmad, warning for madx
        if self.Transport.convprops.gmadoutput:
            self.Transport.machine.AddRCol(name=elementid, length=lenInM, xsize=aperx_in_metres, ysize=apery_in_metres,
                                           material=collimatorMaterial)
        elif self.Transport.convprops.madxoutput:
            self.Writer.DebugPrintout('\tWarning, MadX Builder does not have RCOL')

        debugstring = '\tCollimator, x aperture = ' + _np.str(aperx_in_metres) \
                      + ' m, y aperture = ' + _np.str(apery_in_metres) + ' m.'
        self.Writer.DebugPrintout(debugstring)
        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Collimator ' + elementid + ', length= ' + _np.str(lenInM)\
                      + ' m, xsize= ' + _np.str(_np.round(aperx_in_metres, 4))
        debugstring += ' m, ysize= ' + _np.str(_np.round(apery_in_metres, 4)) + ' m.'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Acceleration(self, linedict):
        """
        A Function that writes the properties of an acceleration element
        Only RF added for gmad, not for madx!
        """
        # Redundant function until comments and /or acceleration components can be handled

        accdata = linedict['data']
        acclen = linedict['length']
        e_gain = linedict['voltage']

        # TODO add phase_lag and wavelength to BDSIM

        # If zero length then start of a sequence, save total accelerating voltage
        if acclen == 0.0:
            self.Transport.convprops.isAccSequence = True
            self.Transport.machineprops._totalAccVoltage = e_gain
            self.Transport.machineprops._e_gain_prev = 0.0  # start at 0
            return

        if self.Transport.convprops.isAccSequence:
            # voltage means voltage relative to the end of this segment
            e_rel_gain = e_gain - self.Transport.machineprops._e_gain_prev
            self.Transport.machineprops._e_gain_prev = e_gain  # store for next segment
            if e_gain == 1.0:  # end of sequence
                self.Transport.convprops.isAccSequence = False

            e_gain = e_rel_gain * self.Transport.machineprops._totalAccVoltage

        gradient = e_gain * (self.Transport.scale[self.Transport.units['p_egain'][0]] / 1e6)
        gradient /= (acclen * self.Transport.scale[self.Transport.units['element_length'][0]])  # gradient in MV/m

        self.Transport.machineprops.rf += 1
        elname = "ACC" + _np.str(self.Transport.machineprops.rf)

        # only call for gmad, warning for madx
        if self.Transport.convprops.gmadoutput:
            self.Transport.machine.AddRFCavity(name=elname, length=acclen, gradient=gradient)
        elif self.Transport.convprops.madxoutput:
            self.Writer.DebugPrintout('\tWarning, MadX Builder does not have RF Cavity')

        # Update beam parameters
        self.Transport = _General.UpdateMomentumFromEnergy(self.Transport, self.Transport.beamprops.k_energy + e_gain)

        # Commented out untested code
        # if len(accdata) == 2:  # Newer case with multiple elements
        # self._acc_sequence(line)
        if len(accdata) == 4:     # Older case for single element
            phase_lag = linedict['phase_lag']
            wavel = linedict['wavel']

            # Write to file
            accline = '! An accelerator element goes here of length ' + \
                      _np.str(acclen) + ' ' + self.Transport.units['element_length'] + ', \n'
            accline2 = '! with an energy gain of ' + _np.str(e_gain) + ' ' + \
                       self.Transport.units['p_egain'] + ', phase lag of ' + _np.str(phase_lag) + ' degrees, \n'
            accline3 = '! and a wavelength of ' + _np.str(wavel) + ' ' + self.Transport.units['bunch_length'] + '. \n'

    def Sextupole(self, linedict):
        sextudata = linedict['data']
        length = sextudata[0]        # First three non-blanks must be the entries in a specific order.
        field_at_tip = sextudata[1]  # Field in TRANSPORT units
        pipe_rad = sextudata[2]      # Pipe Radius In TRANSPORT units

        field_in_Gauss = field_at_tip * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  # Convert to Tesla

        pipe_in_metres = pipe_rad * _General.ScaleToMeters(self.Transport, 'bend_vert_gap')
        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        field_gradient = (2*field_in_Tesla / pipe_in_metres**2) / self.Transport.beamprops.brho  # K2 in correct units

        self.Transport.machineprops.sextus += 1
        elementid = ''
        if self.Transport.convprops.keepName:
            elementid = linedict['name']
        if not elementid:  # check on empty string
            elementid = 'SEXT'+_np.str(self.Transport.machineprops.sextus)

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddSextupole(name=elementid, length=lenInM, k2=_np.round(field_gradient, 4))

        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Sextupole ' + elementid + ', length ' + _np.str(lenInM) + \
                      ' m, k2 ' + _np.str(_np.round(field_gradient, 4)) + ' T/m^2'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Solenoid(self, linedict):
        soledata = linedict['data']
        length = soledata[0]  # First three non-blanks must be the entries in a specific order.
        field = soledata[1]   # Field in TRANSPORT units

        field_in_Gauss = field * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  # Convert to Tesla

        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        self.Transport.machineprops.solenoids += 1
        elementid = ''
        if self.Transport.convprops.keepName:
            elementid = linedict['name']
        if not elementid:  # check on empty string
            elementid = 'SOLE'+_np.str(self.Transport.machineprops.solenoids)

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddSolenoid(name=elementid, length=lenInM, ks=_np.round(field_in_Tesla, 4))

        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Solenoid ' + elementid + ', length ' + _np.str(lenInM) + \
                      ' m, ks ' + _np.str(_np.round(field_in_Tesla, 4)) + ' T'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Printline(self, linedict):
        number = linedict['data'][0]
        self.Writer.DebugPrintout('\tTRANSPORT control line,')
        try:
            number = _np.float(number)
            if number == 48:
                self.Transport.machineprops.benddef = False
                self.Writer.DebugPrintout('\t48. Switched Dipoles to Angle definition.')
            elif number == 47:
                self.Transport.machineprops.benddef = True
                self.Writer.DebugPrintout('\t47. Switched Dipoles to field definition.')
            elif number == 19:
                if _General.CheckSingleLineOutputApplied(self.Transport.convprops.filename):
                    self.Transport.convprops.singleLineOptics = True
                self.Writer.DebugPrintout('\t19. Optics output switched to single line per element.')
            else:
                self.Writer.DebugPrintout('\tCode 13. ' + _np.str(number) + ' handling not implemented.')
        except ValueError:
            pass

    def Correction(self, linedict):
        if self.Transport.convprops.correctedbeamdef:
            self.Writer.DebugPrintout('\tNot Correction to original beam definition')
            return
        # Check if the previous line was the original beam definition and not an rms update
        if linedict['prevlinenum'] == 1.0 and not linedict['isAddition'] and self.Transport.convprops.beamdefined:
            self.Transport.convprops.correctedbeamdef = True

        correctiondata = linedict['data']
        if len(correctiondata) >= 15:  # 15 sigma elements
            sigma21 = correctiondata[0]
            sigma43 = correctiondata[5]
        else:
            self.Writer.DebugPrintout('\tLength of correction line is incorrect')
            return

        emittoverbeta = self.Transport.beamprops.SigmaXP**2 * (1 - sigma21**2)
        emittbeta = self.Transport.beamprops.SigmaX**2
        betx = _np.sqrt(emittbeta / emittoverbeta)
        emitx = emittbeta / betx
        slope = sigma21 * self.Transport.beamprops.SigmaXP / self.Transport.beamprops.SigmaX
        alfx = -1.0 * slope * betx

        self.Transport.beamprops.betx = betx
        self.Transport.beamprops.emitx = emitx / 1000.0
        self.Transport.beamprops.alfx = alfx

        emittoverbeta = self.Transport.beamprops.SigmaYP**2 * (1 - sigma43**2)
        emittbeta = self.Transport.beamprops.SigmaY**2
        bety = _np.sqrt(emittbeta / emittoverbeta)
        emity = emittbeta / bety
        slope = sigma43 * self.Transport.beamprops.SigmaYP / self.Transport.beamprops.SigmaY
        alfy = -1.0 * slope * bety

        self.Transport.beamprops.bety = bety
        self.Transport.beamprops.emity = emity / 1000.0
        self.Transport.beamprops.alfy = alfy

        self.Transport.beamprops.distrType = 'gausstwiss'

        self.Writer.DebugPrintout('\tConverted to:')
        self.Writer.DebugPrintout('\t Beam Correction. Sigma21 = ' + _np.str(sigma21) + ', Sigma43 = ' + _np.str(sigma43) + '.')
        self.Writer.DebugPrintout('\t Beam distribution type now switched to "gausstwiss":')
        self.Writer.DebugPrintout('\t Twiss Params:')
        self.Writer.DebugPrintout('\t BetaX = ' + _np.str(self.Transport.beamprops.betx) + ' ' + self.Transport.units['beta_func'])
        self.Writer.DebugPrintout('\t BetaY = ' + _np.str(self.Transport.beamprops.bety) + ' ' + self.Transport.units['beta_func'])
        self.Writer.DebugPrintout('\t AlphaX = ' + _np.str(self.Transport.beamprops.alfx))
        self.Writer.DebugPrintout('\t AlphaY = ' + _np.str(self.Transport.beamprops.alfy))
        self.Writer.DebugPrintout('\t Emittx = ' + _np.str(self.Transport.beamprops.emitx) + ' ' + self.Transport.units['emittance'])
        self.Writer.DebugPrintout('\t EmittY = ' + _np.str(self.Transport.beamprops.emity) + ' ' + self.Transport.units['emittance'])

    def SpecialInput(self, linedict):
        specialdata = linedict['data']
        self.Writer.DebugPrintout('\tSpecial Input line:')

        if specialdata[0] == 5.0:  # beampiperadius (technically only vertical, but will apply a circle for now)
            self.Writer.DebugPrintout('\tType 5: Vertical half aperture,')
            self.Transport.machineprops.dipoleVertAper = specialdata[1]
            self.Writer.DebugPrintout('\tHalf aperture set to ' + _np.str(specialdata[1]) + '.')
            if self.Transport.machineprops.fringeIntegral == 0:
                self.Transport.machineprops.fringeIntegral = 0.5  # default if a vertical aperture is specified.
                self.Writer.DebugPrintout('FINT/X not set, setting FINT/X to default of 0.5.')
        elif specialdata[0] == 7.0:  # Fringe Field integral
            self.Transport.machineprops.fringeIntegral = specialdata[1]
            self.Writer.DebugPrintout('\tType 7: K1 Fringe field integral,')
            self.Writer.DebugPrintout('\tIntegral set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 14.0:  # Definition of element type code 6.
            if self.Transport.convprops.typeCode6IsTransUpdate:
                self.Transport.convprops.typeCode6IsTransUpdate = False
                typeCode6def = 'Collimator'
            else:
                self.Transport.convprops.typeCode6IsTransUpdate = True
                typeCode6def = 'Transform Update'
            self.Writer.DebugPrintout('\tType 14: Type code 6 definition,')
            self.Writer.DebugPrintout('\tDefinition set to ' + typeCode6def + '.')
        elif specialdata[0] == 16.0:  # X0 offset
            self.Transport.beamprops.X0 = specialdata[1]
            self.Writer.DebugPrintout('\tType 16: X0 beam offset,')
            self.Writer.DebugPrintout('\tOffset set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 17.0:  # Y0 offset
            self.Transport.beamprops.Y0 = specialdata[1]
            self.Writer.DebugPrintout('\tType 17: Y0 beam offset,')
            self.Writer.DebugPrintout('\tOffset set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 18.0:  # Z0 offset
            self.Transport.beamprops.Z0 = specialdata[1]
            self.Writer.DebugPrintout('\tType 18: Z0 beam offset,')
            self.Writer.DebugPrintout('\tOffset set to ' + _np.str(specialdata[1]) + '.')
        else:
            self.Writer.DebugPrintout('\tCode type not yet supported, or unknown code type.')

    def UnitChange(self, linedict):
        """
        Function to change the units (scaling) of various parameters.
        """
        label = linedict['label']
        number = linedict['number']
        if label == 'CM' or label == 'MM' or label == 'UM' or label == 'NM':
            label = label.lower()
        # Convert Energy Unit Cases:
        if label == 'EV':
            label = 'eV'
        if label == 'KEV':
            label = 'keV'
        if label == 'MEV':
            label = 'MeV'
        if label == 'GEV':
            label = 'GeV'
        if label == 'TEV':
            label = 'TeV'

        debugstring2 = '\tConverted to ' + label

        if _np.float(number) == 1:  # Horizontal and vertical beam size
            self.Transport.units['x'] = label
            self.Transport.units['y'] = label
            self.Transport.units['bend_vert_gap'] = label
            # self.units['pipe_rad'] = label
            debugstring1 = '\tType 1: Horizontal and vertical beam extents, and magnet apertures,'

        elif _np.float(number) == 2:  # Horizontal and vertical divergence
            self.Transport.units['xp'] = label
            self.Transport.units['yp'] = label
            debugstring1 = '\tType 2: Horizontal and vertical angles,'

        elif _np.float(number) == 3:  # Bending Magnet Gap
            self.Transport.units['y'] = label
            self.Transport.units['bend_vert_gap'] = label
            debugstring1 = '\tType 3: Vertical (only) beam extent and magnet aperture,'

        elif _np.float(number) == 4:  # Vertical Divergence ONLY
            self.Transport.units['yp'] = label
            debugstring1 = '\tType 4: Vertical (only) beam angle,'

        elif _np.float(number) == 5:  # Pulsed Beam Length
            self.Transport.units['bunch_length'] = label
            debugstring1 = '\tType 5: Bunch length,'

        elif _np.float(number) == 6:  # Momentum Spread
            self.Transport.units['momentum_spread'] = label  # Percent
            debugstring1 = '\tType 6: Momentum spread,'

        elif _np.float(number) == 7:  # Bend/pole face rotation
            debugstring1 = '\tType 7: Bend and poleface rotation angles,'
            debugstring2 = '\tCONVERTION NOT IMPLEMENTED YET.'
            pass

        elif _np.float(number) == 8:  # Element Length
            self.Transport.units['element_length'] = label
            debugstring1 = '\tType 8: Element length,'

        elif _np.float(number) == 9:  # Magnetic Field
            self.Transport.units['magnetic_fields'] = label
            debugstring1 = '\tType 9: Magnetic Fields,'

        elif _np.float(number) == 10:  # Mass
            debugstring1 = '\tType 10: Mass,'
            debugstring2 = '\tCONVERTION NOT IMPLEMENTED YET.'
            pass

        elif _np.float(number) == 11:  # Momentum / energy gain during acc.
            self.Transport.units['p_egain'] = label
            debugstring1 = '\tType 11: Momentum and accelerator energy gain,'
        else:
            # default output
            debugstring1 = '\tCode type not yet supported, or unknown code type.'
            debugstring2 = ''

        self.Writer.DebugPrintout('\tUnit change line:')
        self.Writer.DebugPrintout(debugstring1)
        self.Writer.DebugPrintout(debugstring2)

    def TransformUpdate(self, linedict):
        if linedict['elementnum'] == 6.0:
            errorline = '\tElement is either a transform update or a collimator. The type code 6 definition'
            errorline2 = '\thas not been switched to collimators, therefore nothing will be done for this element.'
            self.Writer.DebugPrintout(errorline)
            self.Writer.DebugPrintout(errorline2)
