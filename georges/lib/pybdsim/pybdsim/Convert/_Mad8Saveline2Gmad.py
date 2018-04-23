from .. import Builder as _Builder
from .. import Beam as _Beam
from  pymad8.Output import Saveline as _Saveline

# pybdsim.Convert.Mad8Saveline2Gmad('../ebds.saveline', 'ilc.gmad', ignore_zero_length_items=False)

def Mad8Saveline2Gmad(input, output_file_name, start_name=None, end_name=None, ignore_zero_length_items=True,
                      samplers='all', aperture_dict={}, collimator_dict="collimators.dat", beam_pipe_radius=0.2,
                      verbose=False, beam=True, optics=True, loss=True):
    items_omitted = []  # Store any items that have been skipped over.
    fake_length = 1e-6

    if type(input) == str:
        mad8 = _Saveline.Loader(input)
    else:
        mad8 = input

    if type(collimator_dict) == str :
        collimator_obj = Mad8SavelineCollimatorDatabase(collimator_dict)
        collimator_dict = collimator_obj.GetDict()
        xsize_dict = 'xsize'
        ysize_dict = 'ysize'
    else:  # Necessary as the dictionaries use different cases.
        xsize_dict = 'XSIZE'
        ysize_dict = 'YSIZE'
    
    ilc = _Builder.Machine()

    for element in mad8.elementList:
        element_type = mad8.elementDict[element].keys()[0]
        element_properties = mad8.elementDict[element].values()[0]
        length = element_properties['L'] if 'L' in element_properties else fake_length  # Everything has a length.
        zero_length = True if length < 1e-9 else False

        if element_type == 'LINE':
            for name in element_properties:
                sequence_element_type = mad8.elementDict[name].keys()[0]
                sequence_element_properties = mad8.elementDict[name].values()[0]
                # Everything has a length.
                length = sequence_element_properties['L'] if 'L' in sequence_element_properties else fake_length

                construct_sequence(ilc, name, sequence_element_type, sequence_element_properties, length,
                                   ignore_zero_length_items, items_omitted, zero_length, aperture_dict, collimator_dict,
                                   beam_pipe_radius, xsize_dict, ysize_dict, verbose)

    ilc.AddSampler(samplers)

    if beam:  # Use ILC defaults
        beam = _Beam(particletype='e-', energy=250, distrtype='gausstwiss')

        beam._SetBetaX(71.482996720846, unitsstring='m')
        beam._SetBetaY(39.603963961143, unitsstring='m')
        beam._SetAlphaX(-1.562625067837)
        beam._SetAlphaY(1.283201237991)
        beam.SetZ0(-1, unitsstring='um')
        beam.SetT0(1e-15)
        beam._SetEmittanceX(emitx=2.044e-11,unitsstring='m')
        beam._SetEmittanceY(emity=8.176e-14, unitsstring='m')

        ilc.AddBeam(beam)

    if optics:
        ilc.options.SetPhysicsList('QGSP_BERT')

    if loss:
        ilc.options.SetStopTracks(False)
            
    ilc.options.SetNGenerate(1000)
    ilc.Write(output_file_name)


def construct_sequence(ilc, element, element_type, element_properties, length, ignore_zero_length_items,
                       items_omitted, zero_length, aperture_dict, collimator_dict, beam_pipe_radius, xsize_dict, ysize_dict, verbose):
    kws = {}
    if element_type == 'DRIFT':
            ilc.AddDrift(element, length, **kws)

    elif element_type == 'VKICKER':
        kws['kick'] = element_properties['KICK'] if 'KICK' in element_properties else 0

        ilc.AddVKicker(element, length, **kws)

    elif element_type == 'HKICKER':
        kws['kick'] = element_properties['KICK'] if 'KICK' in element_properties else 0

        ilc.AddHKicker(element, length, **kws)

    elif element_type == 'MARKER':
        if ignore_zero_length_items:
            items_omitted.append(element)
        else:
            ilc.AddMarker(element)

    elif element_type == 'SBEND':
        angle = element_properties['ANGLE'] if 'ANGLE' in element_properties else 0
        k1 = (element_properties['K1'] / length) if 'K1' in element_properties else 0

        kws['hgap'] = element_properties['HGAP'] if 'HGAP' in element_properties else 0
        kws['fintx'] = element_properties['FINTX'] if 'FINTX' in element_properties else 0
        kws['tilt'] = element_properties['TILT'] if 'TILT' in element_properties else 0
        kws['fint'] = element_properties['FINT'] if 'FINT' in element_properties else 0
        kws['e1'] = element_properties['E1'] if 'E1' in element_properties else 0
        kws['e2'] = element_properties['E2'] if 'E2' in element_properties else 0

        ilc.AddDipole(element, 'sbend', length, angle, k1=k1, **kws)

    elif element_type == 'RBEND':
        angle = element_properties['ANGLE'] if 'ANGLE' in element_properties else 0

        ilc.AddDipole(element, 'rbend', angle, **kws)

    elif element_type == 'QUADRUPOLE':
        k1 = element_properties['K1'] if 'K1' in element_properties else 0

        #kws['aperture'] = element_properties['APERTURE'] if 'APERTURE' in element_properties else 0

        ilc.AddQuadrupole(element, length, k1, **kws)

    elif element_type == 'SEXTUPOLE':
        k2 = element_properties['K2'] if 'K2' in element_properties else 0

        kws['tilt'] = element_properties['TILT'] if 'TILT' in element_properties else 0

        #  Aperture not currently supported
        #kws['aperture'] = element_properties['APERTURE'] if 'APERTURE' in element_properties else 0 aperX / aperY

        #ilc.AddSextupole(element, length, k2, **kws)
        ilc.AddDrift(element,length)
    elif element_type == 'OCTUPOLE':
        k3 = element_properties['K3'] if 'K3' in element_properties else 0

        #  Aperture not currently supported.
        #kws['aperture'] = element_properties['APERTURE'] if 'APERTURE' in element_properties else 0

        # ilc.AddOctupole(element, length, k3, **kws)
        ilc.AddDrift(element, length)

    elif element_type == 'MULTIPOLE':  # Not fully implemented, but ready to go when AddMultipole works
        tilt = element_properties['TILT'] if 'TILT' in element_properties else 0
        knl = ()
        kns = ()

        for i in range(1, 7):  # Range is fixed to between 1 - 6 as far as I can tell, might be more though
            key_l = 'K%sL' % i
            key_s = 'K%sS' % i
            temp = element_properties[key_l] if key_l in element_properties else 0
            temp2 = element_properties[key_s] if key_s in element_properties else 0
            knl = knl + (temp,)
            kns = kns + (temp2,)

        kws['lrad'] = element_properties['LRAD'] if 'LRAD' in element_properties else 0
        #  Aperture not currently supported
        #kws['aperture'] = element_properties['APERTURE'] if 'APERTURE' in element_properties else 0

        # ilc.AddMultipole(element, length, knl, ksl=kns, tilt=tilt, **kws)  # ksl should be kns to keep format
        ilc.AddMarker(element)

    elif element_type == 'ECOLLIMATOR': # Need to check logic precedence here with Boog.
        if element in collimator_dict:
            kws['material'] = collimator_dict[element]['bdsim_material'] if 'bdsim_material' in collimator_dict else 'Copper'
            xsize = collimator_dict[element][xsize_dict] if xsize_dict in collimator_dict[element] else beam_pipe_radius
            ysize = collimator_dict[element][ysize_dict] if ysize_dict in collimator_dict[element] else beam_pipe_radius
#            angle = collimator_dict[element]['ANGLE'] if 'ANGLE' in collimator_dict[element] else 0.0
        else:
            kws['material'] = 'Copper'
            xsize = element_properties[xsize_dict] if xsize_dict in element_properties else beam_pipe_radius
            ysize = element_properties[ysize_dict] if ysize_dict in element_properties else beam_pipe_radius
#            angle = element_properties['ANGLE'] if 'ANGLE' in element_properties else 0

        ilc.AddECol(element, length, xsize, ysize, **kws)

    elif element_type == 'RCOLLIMATOR':
        if element in collimator_dict:
            kws['material'] = collimator_dict[element]['bdsim_material'] if 'bdsim_material' in collimator_dict else 'Copper'
            xsize = collimator_dict[element][xsize_dict] if xsize_dict in collimator_dict[element] else beam_pipe_radius
            ysize = collimator_dict[element][ysize_dict] if ysize_dict in collimator_dict[element] else beam_pipe_radius
 #           angle = collimator_dict[element]['ANGLE'] if 'ANGLE' in collimator_dict[element] else 0.0
        else:
            kws['material'] = 'Copper'
            xsize = element_properties[xsize_dict] if xsize_dict in element_properties else beam_pipe_radius
            ysize = element_properties[ysize_dict] if ysize_dict in element_properties else beam_pipe_radius
#            angle = element_properties['ANGLE'] if 'ANGLE' in element_properties else 0

        ilc.AddRCol(element, length, xsize, ysize, **kws)

    elif element_type == 'WIRE':
        if ignore_zero_length_items and zero_length:
            items_omitted.append(element)
        elif (not ignore_zero_length_items) and zero_length:
            ilc.AddMarker(element)
            if verbose:
                print element, ' -> marker instead of wire.'
        else:
            ilc.AddDrift(element, length, **kws)

    elif element_type == 'INSTRUMENT':
        if ignore_zero_length_items and zero_length:
            items_omitted.append(element)
        elif (not ignore_zero_length_items) and zero_length:
            ilc.AddMarker(element)
            if verbose:
                print element, ' -> marker instead of instrument.'
        else:
            ilc.AddDrift(element, length, **kws)

    elif element_type == 'MONITOR':
        if ignore_zero_length_items and zero_length:
            items_omitted.append(element)
        elif (not ignore_zero_length_items) and zero_length:
            ilc.AddMarker(element)
            if verbose:
                print element, ' -> marker instead of monitor.'
        else:
            ilc.AddDrift(element, length, **kws)

    elif element_type == 'LCAVITY':
        deltae = element_properties['DELTAE'] * 1000 if 'DELTAE' in element_properties else 0  # MeV
        gradient = deltae / length

        #kws['aperture'] = element_properties['APERTURE'] if 'APERTURE' in element_properties else 0
        kws['freq'] = element_properties['FREQ'] if 'FREQ' in element_properties else 0
        kws['phi0'] = element_properties['PHI0'] if 'PHI0' in element_properties else 0
        kws['eloss'] = element_properties['ELOSS'] if 'ELOSS' in element_properties else 0

        ilc.AddRFCavity(element, length, gradient, **kws)


class Mad8SavelineCollimatorDatabase: 
    '''
    Load collimator file into memory and functions to open and manipulate collimator system
    c = Mad8CollimatorDataBase(fileName)
    fileName = "collimator.dat"
    file format
    <name> <angle> <length> <x_size/2> <ysize/2> <material>
    <length> includes the tapers, so entire length along machine 
    '''

    def __init__(self,collimatorFileName) :
        self.collimatorFileName = collimatorFileName
        self.LoadCollimatorDb(self.collimatorFileName) 
        
    def LoadCollimatorDb(self,collimatorFileName) : 
        f = open(collimatorFileName) 

        inx = 0 

        self._coll = {}
        for l in f : 
            t = l.split()
            name     = t[0]
            angle    = float(t[1])
            length   = float(t[2])
            xsize    = float(t[3])
            ysize    = float(t[4]) 
            material = t[5]
            inx = inx + 1 

            d = {'index':inx, 'angle':angle, 'l':length, 'xsize':xsize,
                 'ysize':ysize, 'bdsim_material':material}

            self._coll[name] = d     
        
    def OpenCollimators(self,openHalfSizeX=0.1, openHalfSizeY=0.1) : 
        for c in self._coll.keys() :
            self._coll[c]['xsize'] = openHalfSizeX
            self._coll[c]['ysize'] = openHalfSizeY
    
    def GetCollimators(self) : 
        return self._coll.keys()

    def GetCollimator(self, name) : 
        return self._coll[name]

    def GetDict(self) : 
        return self._coll
