from copy import deepcopy as _deepcopy
import numpy as _np
import re as _re
import pymadx as _pymadx
import warnings as _warnings
from .. import Builder as _Builder
from ..Options import Options as _Options
from .. import Beam as _Beam
from .. import _General
from .. import XSecBias

_requiredKeys = frozenset([
    'L', 'ANGLE', 'KSI', 'K1L', 'K2L', 'K3L', 'K4L', 'K5L',
    'K1SL', 'K2SL', 'K3SL', 'K4SL', 'K5SL', 'K6SL',
    'TILT', 'KEYWORD', 'ALFX', 'ALFY', 'BETX', 'BETY',
    'VKICK', 'HKICK', 'E1', 'E2', 'FINT', 'FINTX'])

_lFake = 1e-6 # fake length for thin magnets

def ZeroMissingRequiredColumns(tfsinstance):
    """
    Sets any missing required columns to zero.  Warns user when doing so.
    """
    missingColumns = [key for key
                      in _requiredKeys if key not in tfsinstance.columns]

    if not missingColumns:
        return

    for column in missingColumns:
        tfsinstance.columns.append(column)
        for key, data in tfsinstance.data.iteritems():
            data.append(0.0)

    missingColsString = ", ".join(["\"{}\"".format(col)
                                   for col in missingColumns])
    msg = ("Columns missing from TFS: {}.  All have been set"
           " to zero.").format(missingColsString)
    print(msg)


def MadxTfs2Gmad(tfs, outputfilename, startname=None, stopname=None, stepsize=1,
                 ignorezerolengthitems=True, thinmultipoles=True, samplers='all',
                 aperturedict={},
                 collimatordict={},
                 userdict={},
                 verbose=False, beam=True, flipmagnets=None, usemadxaperture=False,
                 defaultAperture='circular',
                 biases=None,
                 allelementdict={},
                 optionsDict={},
                 linear=False,
                 overwrite=True,
                 allNamesUnique=False,
                 gaussmatrix=False):
    """
    **MadxTfs2Gmad** convert a madx twiss output file (.tfs) into a gmad tfs file for bdsim

    +-------------------------------+-------------------------------------------------------------------+
    | **tfs**                       | path to the input tfs file or pymadx.Data.Tfs instance            |
    +-------------------------------+-------------------------------------------------------------------+
    | **outputfilename**            | requested output file                                             |
    +-------------------------------+-------------------------------------------------------------------+
    | **startname**                 | the name (exact string match) of the lattice element to start the |
    |                               | machine at this can also be an integer index of the element       |
    |                               | sequence number in madx tfs.                                      |
    +-------------------------------+-------------------------------------------------------------------+
    | **stopname**                  | the name (exact string match) of the lattice element to stop the  |
    |                               | machine at this can also be an integer index of the element       |
    |                               | sequence number in madx tfs.                                      |
    +-------------------------------+-------------------------------------------------------------------+
    | **stepsize**                  | the slice step size. Default is 1, but -1 also useful for         |
    |                               | reversed line.                                                    |
    +-------------------------------+-------------------------------------------------------------------+
    | **ignorezerolengthitems**     | nothing can be zero length in bdsim as real objects of course     |
    |                               | have some finite size.  Markers, etc are acceptable but for large |
    |                               | lattices this can slow things down. True allows to ignore these   |
    |                               | altogether, which doesn't affect the length of the machine.       |
    +-------------------------------+-------------------------------------------------------------------+
    | **thinmultipoles**            | will convert thin multipoles to ~1um thick finite length          |
    |                               | multipoles with upscaled k values - experimental feature          |
    +-------------------------------+-------------------------------------------------------------------+
    | **samplers**                  | can specify where to set samplers - options are None, 'all', or a |
    |                               | list of names of elements (normal python list of strings). Note   |
    |                               | default 'all' will generate separate outputfilename_samplers.gmad |
    |                               | with all the samplers which will be included in the main .gmad    |
    |                               | file - you can comment out the include to therefore exclude all   |
    |                               | samplers and retain the samplers file.                            |
    +-------------------------------+-------------------------------------------------------------------+
    | **aperturedict**              | Aperture information. Can either be a dictionary of dictionaries  |
    |                               | with the the first key the exact name of the element and the      |
    |                               | daughter dictionary containing the relevant bdsim parameters as   |
    |                               | keys (must be valid bdsim syntax). Alternatively, this can be a   |
    |                               | pymadx.Aperture instance that will be queried.                    |
    +-------------------------------+-------------------------------------------------------------------+
    | **collimatordict**            | A dictionary of dictionaries with collimator information keys     |
    |                               | should be exact string match of element name in tfs file value    |
    |                               | should be dictionary with the following keys:                     |
    |                               | "bdsim_material"   - the material                                 |
    |                               | "angle"            - rotation angle of collimator in radians      |
    |                               | "xsize"            - x full width in metres                       |
    |                               | "ysize"            - y full width in metres                       |
    +-------------------------------+-------------------------------------------------------------------+
    | **userdict**                  | A python dictionary the user can supply with any additional       |
    |                               | information for that particular element. The dictionary should    |
    |                               | have keys matching the exact element name in the Tfs file and     |
    |                               | contain a dictionary itself with key, value pairs of parameters   |
    |                               | and values to be added to that particular element.                |
    +-------------------------------+-------------------------------------------------------------------+
    | **verbose**                   | Print out lots of information when building the model.            |
    +-------------------------------+-------------------------------------------------------------------+
    | **beam**                      | True \| False - generate an input gauss Twiss beam based on the   |
    |                               | values of the twiss parameters at the beginning of the lattice    |
    |                               | (startname) NOTE - we thoroughly recommend checking these         |
    |                               | parameters and this functionality is only for partial convenience |
    |                               | to have a model that works straight away.                         |
    +-------------------------------+-------------------------------------------------------------------+
    | **flipmagnets**               | True \| False - flip the sign of all k values for magnets - MADX  |
    |                               | currently tracks particles agnostic of the particle charge -      |
    |                               | BDISM however, follows their manual definition strictly -         |
    |                               | positive k -> horizontal focussing for positive particles         |
    |                               | therefore, positive k -> vertical focussing for negative          |
    |                               | particles. Use this flag to flip the sign of all magnets.         |
    +-------------------------------+-------------------------------------------------------------------+
    | **usemadxaperture**           | True \| False - use the aperture information in the TFS file if   |
    |                               | APER_1 and APER_2 columns exist.  Will only set if they're        |
    |                               | non-zero.                                                         |
    +-------------------------------+-------------------------------------------------------------------+
    | **defaultAperture**           | The default aperture model to assume if none is specified.        |
    +-------------------------------+-------------------------------------------------------------------+
    | **biases**                    | Optional list of bias objects to be defined in own _bias.gmad     |
    |                               | file.  These can then be attached either with allelementdict for  |
    |                               | all components or userdict for individual ones.                   |
    +-------------------------------+-------------------------------------------------------------------+
    | **allelementdict**            | Dictionary of parameter/value pairs to be written to all          |
    |                               | components.                                                       |
    +-------------------------------+-------------------------------------------------------------------+
    | **optionsDict**               | Optional dictionary of general options to be written to the       |
    |                               | bdsim model options.                                              |
    +-------------------------------+-------------------------------------------------------------------+
    | **linear**                    | Only linear optical components                                    |
    +-------------------------------+-------------------------------------------------------------------+
    | **overwrite**                 | Do not append an integer to the base file name if it already      |
    |                               | exists.  Instead overwrite the files.                             |
    +-------------------------------+-------------------------------------------------------------------+
    | **allNamesUnique              | Treat every row in the TFS file/instance as a unique element.     |
    |                               | This makes it easier to edit individual components as they are    |
    |                               | guaranteed to appear only once in the entire resulting GMAD       |
    |                               | lattice.                                                          |
    +-------------------------------+-------------------------------------------------------------------+

    Example:

    >>> a,o = pybdsim.Convert.MadxTfs2Gmad('twiss.tfs', 'mymachine')

    In normal mode:
    Returns Machine, [omittedItems]

    In verbose mode:
    Returns Machine, Machine, [omittedItems]

    Returns two pybdsim.Builder.Machine instances. The first desired full conversion.  The second is
    the raw conversion that's not split by aperture. Thirdly, a list of the names of the omitted items
    is returned.
    """

    # machine instance that will be added to
    a = _Builder.Machine() # raw converted machine
    b = _Builder.Machine() # final machine, split with aperture

    # test whether filepath or tfs instance supplied
    madx = _pymadx.Data.CheckItsTfs(tfs)

    factor = 1
    if flipmagnets != None:
        factor = -1 if flipmagnets else 1  #flipping magnets
    else:
        if madx.header.has_key('PARTICLE'):
            particleName = madx.header['PARTICLE']
            if particleName == "ELECTRON":
                flipmagnets = True
                print('Detected electron in TFS file - changing flipmagnets to True')


    # If we have collimators but no collimator dict then inform that
    # they will be converted to drifts.  should really check
    # tfs[startname..]
    if "APERTYPE" in madx.columns:
        if (("RCOLLIMATOR" in madx.GetColumn("APERTYPE")
            or "ECOLLIMATOR" in madx.GetColumn("APERTYPE"))
            and not collimatordict):
            warning.warn("No collimatordict provided.  ALL collimators"
                         " will be converted to DRIFTs.")

    if isinstance(biases, XSecBias.XSecBias):
        a.AddBias(biases)
        b.AddBias(biases)
    elif biases is not None:
        try:
            [a.AddBias(bias) for bias in biases]
            [b.AddBias(bias) for bias in biases]
        except TypeError:
            raise TypeError("Unknown biases!  Biases must be a XSecBias"
                            "instance or an iterable of XSecBias instances.")

    # define utility function that does conversion
    def AddSingleElement(item, a, aperModel=None):
        # a is a pybdsim.Builder.Machine instance
        # if it's already a prepared element, just append it
        if type(item) == _Builder.Element:
            a.Append(item)
            return

        kws = {} # ensure empty
        # deep copy as otherwise allelementdict gets irreperably changed!
        kws = _deepcopy(allelementdict)
        if verbose:
            print('starting key word arguments from all element dict')
            print(kws)
        
        if aperModel != None:
            kws.update(aperModel)

        name  = item['NAME']
        # remove special characters like $, % etc 'reduced' name - rname:
        rname = pybdsim._General.PrepareReducedName(name
                                                    if not allNamesUnique
                                                    else item["UNIQUENAME"])
        t     = item['KEYWORD']
        l     = item['L']
        ang   = item['ANGLE']
        tilt  = item['TILT']

        if tilt != 0:
            kws['tilt'] = tilt

        # append any user defined parameters for this element into the
        # kws dictionary

        # name appears in the madx.  try this first.
        if name in userdict:
            kws.update(userdict[name])
        elif rname in userdict: # rname appears in the gmad
            kws.update(userdict[rname])

        if verbose:
            print('Full set of key word arguments:')
            print(kws)

        if t == 'DRIFT':
            a.AddDrift(rname,l,**kws)
        elif t == 'HKICKER':
            if verbose:
                print('HICKER',rname)
            hkick = item['HKICK'] * factor
            if not zerolength:
                kws['l'] = l
            a.AddHKicker(rname,hkick=hkick,**kws)
        elif t == 'VKICKER':
            if verbose:
                print('VKICKER',rname)
            vkick = item['VKICK'] * factor
            if not zerolength:
                kws['l'] = l
            a.AddVKicker(rname,vkick=vkick,**kws)
        elif t == 'KICKER':
            if verbose:
                print('KICKER',rname)
            hkick = item['HKICK'] * factor
            vkick = item['VKICK'] * factor
            if not zerolength:
                kws['l'] = l
            a.AddKicker(rname,hkick=hkick,vkick=vkick,**kws)
        elif t == 'TKICKER':
            if verbose:
                print('TKICKER',rname)
            hkick = item['HKICK'] * factor
            vkick = item['VKICK'] * factor
            if not zerolength:
                kws['l'] = l
            a.AddTKicker(rname,hkick=hkick,vkick=vkick,**kws)
        elif t == 'INSTRUMENT':
            #most 'instruments' are just markers
            if zerolength and not ignorezerolengthitems:
                a.AddMarker(rname)
                if verbose:
                    print(name,' -> marker instead of instrument')
            else:
                a.AddDrift(rname,l,**kws)
        elif t == 'MARKER':
            if not ignorezerolengthitems:
                a.AddMarker(rname)
        elif t == 'MONITOR':
            #most monitors are just markers
            if zerolength and not ignorezerolengthitems:
                a.AddMarker(rname)
                if verbose:
                    print(name,' -> marker instead of monitor')
            else:
                a.AddDrift(rname,l,**kws)
        elif t == 'MULTIPOLE':
            k1  = item['K1L'] * factor
            k2  = item['K2L'] * factor
            k3  = item['K3L'] * factor
            k4  = item['K4L'] * factor
            k5  = item['K5L'] * factor
            k6  = item['K6L'] * factor
            k1s = item['K1SL'] * factor
            k2s = item['K2SL'] * factor
            k3s = item['K3SL'] * factor
            k4s = item['K4SL'] * factor
            k5s = item['K5SL'] * factor
            k6s = item['K6SL'] * factor
            tilt= item['TILT']

            if thinmultipoles and not linear :
                a.AddThinMultipole(rname,
                                   knl=(k1, k2, k3, k4, k5, k6),
                                   ksl=(k1s, k2s, k3s, k4s, k5s, k6s),
                                   **kws)
            elif zerolength and not ignorezerolengthitems:
                a.AddMarker(rname)
                if verbose:
                    print(name,' -> marker instead of multipole')
            else:
                a.AddDrift(rname,l,**kws)
        elif t == 'OCTUPOLE':
            if linear :
                k3 = 0.0
            else :
                k3 = item['K3L'] / l * factor
            a.AddOctupole(rname,l,k3=k3,**kws)
        elif t == 'PLACEHOLDER':
            if zerolength:
                a.AddMarker(rname)
                if verbose:
                    print(name,' -> marker instead of placeholder')
            else:
                a.AddDrift(rname,l,**kws)
        elif t == 'QUADRUPOLE':
            k1 = item['K1L'] / l * factor
            a.AddQuadrupole(rname,l,k1=k1,**kws)
        elif t == 'RBEND':
            angle = item['ANGLE']
            e1    = item['E1']
            e2    = item['E2']
            fint  = item['FINT']
            fintx = item['FINTX']
            hgap = item.get('HGAP', 0.0)
            k1l = item["K1L"]

            if (e1 != 0):
                kws['e1'] = e1
            if (e1 != 0):
                kws['e2'] = e2
            if (fint != 0):
                kws['fint']  = fint
            if (fintx != 0):
                kws['fintx'] = fintx
            if (hgap != 0):
                kws['hgap'] = hgap
            if k1l != 0:
                k1 = k1l / l
                kws['k1'] = k1
            a.AddDipole(rname,'rbend',l,angle=angle,**kws)
        elif t == 'RCOLLIMATOR' or t == 'ECOLLIMATOR':
            #only use xsize as only have half gap
            if name in collimatordict:
                #gets a dictionary then extends kws dict with that dictionary
                colld = collimatordict[name]

                #collimator defined by external geometry file
                if 'geometryFile' in colld:
                    kws['geometryFile'] = colld['geometryFile']
                    k = 'outerDiameter' #key
                    if k not in kws:
                        #not already specified via other dictionaries
                        if k in colld:
                            kws[k] = colld[k]
                        else:
                            kws[k] = 1 #ensure there's a default as not in madx
                    # add a general element
                    a.AddElement(rname,l,**kws)
                else:
                    kws['material'] = colld['material']
                    kws['tilt']     = colld['tilt']
                    xsize           = colld['xsize']
                    ysize           = colld['ysize']
                    if verbose:
                        print('collimator xsize ',xsize)
                    #if xsize > 0.1 or ysize > 0.1:
                    kws['outerDiameter'] = max(xsize,ysize)*2.5
                    if t == 'RCOLLIMATOR':
                        a.AddRCol(rname,l,xsize,ysize,**kws)
                    else:
                        a.AddECol(rname,l,xsize,ysize,**kws)
            # if user provided an incomplete collimatordict.
            elif collimatordict != {}:
                msg = ("{} {} not found in collimatordict."
                       " Will instead convert to a DRIFT!".format(t, name))
                _warnings.warn(msg)
                a.AddDrift(rname,l)
            # if user didn't provide a collimatordict at all.
            else:
                a.AddDrift(rname,l)
        elif t == 'RFCAVITY':
            a.AddDrift(rname,l,**kws)
        elif t == 'SBEND':
            angle = item['ANGLE']
            e1    = item['E1']
            e2    = item['E2']
            fint  = item['FINT']
            fintx = item['FINTX']
            k1l   = item['K1L']
            hgap = item.get('HGAP', 0.0)
            if (e1 != 0):
                kws['e1'] = e1
            if (e2 != 0):
                kws['e2'] = e2
            if k1l != 0:
                k1 = k1l / l
                kws['k1'] = k1
            if (fint != 0):
                kws['fint']  = fint
            if (fintx != 0):
                kws['fintx'] = fintx
            if (hgap != 0):
                kws['hgap'] = hgap
            a.AddDipole(rname,'sbend',l,angle=angle,**kws)
        elif t == 'SEXTUPOLE':
            if linear :
                k2 = 0.0
            else :
                k2 = item['K2L'] / l * factor
            a.AddSextupole(rname,l,k2=k2,**kws)
        elif t == 'SOLENOID':
            #ks = item['KSI'] / l
            #a.AddSolenoid(rname,l,ks=ks
            a.AddDrift(rname,l,**kws)
        elif t == 'TKICKER':
            a.AddDrift(rname,l,**kws)
        else:
            print('unknown element type:', t, 'for element named: ', name)
            if zerolength and not ignorezerolengthitems:
                print('putting marker in instead as its zero length')
                a.AddMarker(rname)
            else:
                print('putting drift in instead as it has a finite length')
                a.AddDrift(rname,l)
    # end of utility conversion function

    # check whether it has all the required columns.
    ZeroMissingRequiredColumns(madx)

    if verbose:
        madx.ReportPopulations()

    # check aperture information if supplied
    useTfsAperture = False
    if type(aperturedict) == _pymadx.Data.Aperture:
        useTfsAperture = True
        if verbose:
            aperturedict.ReportPopulations()
    if verbose:
        print('Using pymadx.Apeture instance? --> ',useTfsAperture)

    # keep list of omitted zero length items
    itemsomitted = []

    # iterate through input file and construct machine
    for item in madx[startname:stopname:stepsize]:
        name = item['NAME']
        t    = item['KEYWORD']
        l    = item['L']
        zerolength = True if item['L'] < 1e-9 else False
        if verbose:
            print('zerolength? ',str(name).ljust(20),str(l).ljust(20),' ->',zerolength)
        if madx.ElementPerturbs(item):
            pass #ie proceed normally
        elif zerolength and ignorezerolengthitems:
            itemsomitted.append(name)
            if verbose:
                print('skipping this item')
            continue # skip this item in the for loop
        
        # now deal with aperture
        if useTfsAperture:
            sMid = (item['S']*2 - item['L'] ) * 0.5
            apermodel = _Builder.PrepareApertureModel(aperturedict.GetApertureAtS(sMid), defaultAperture)
            #apermodel = aperturedict.GetApertureForElementNamed(name)
            #print 'Using aperture instance'
            should,lengths,apers = aperturedict.ShouldSplit(item)
            should = False
            if should:
                if verbose:
                    print('Splitting item based on aperture')
                ls = _np.array(lengths)
                if abs(_np.array(lengths).sum() - l) > 1e-6 or (ls < 0).any():
                    print('OH NO!!!')
                    print(l)
                    print(lengths)
                    print(apers)
                    return
                # we should split this item up
                # add it first to the raw machine
                AddSingleElement(item, a)
                # now get the last element - the one that's just been added
                lastelement = a[-1]
                for splitLength,aper in zip(lengths,apers):
                    # append the right fraction with the appropriate aperture
                    # to the 'b' machine
                    apermodel = _Builder.PrepareApertureModel(aper, defaultAperture)
                    print(apermodel)
                    AddSingleElement(lastelement*(splitLength/l), b, apermodel)
            #else:
                #apermodel = _Builder.PrepareApertureModel(apers[0],defaultAperture)
                #print apermodel
                #print len(b)
            AddSingleElement(item, a, apermodel)
            AddSingleElement(item, b, apermodel)
            #    print len(b)
        elif usemadxaperture and name not in aperturedict:
            print('Using aperture in madx tfs file')
            apermodel = _Builder.PrepareApertureModel(item, defaultAperture)
            AddSingleElement(item, a, apermodel)
            AddSingleElement(item, b, apermodel)
        elif item['NAME'] in aperturedict:
            apermodel = _Builder.PrepareApertureModel(aperturedict[name], defaultAperture)
            AddSingleElement(item, a, apermodel)
            AddSingleElement(item, b, apermodel)
        else:
            AddSingleElement(item, a)
            AddSingleElement(item, b)
    # end of for loop

    #add a single marker at the end of the line
    a.AddMarker('theendoftheline')
    b.AddMarker('theendoftheline')

    if (samplers != None):
        a.AddSampler(samplers)
        b.AddSampler(samplers)

    # Make beam file
    if beam:
        bm = MadxTfs2GmadBeam(madx, startname, verbose, gaussmatrix)
        a.AddBeam(bm)
        b.AddBeam(bm)

    options = _Options()
    if (len(optionsDict) > 0):
        options.update(optionsDict) # expand with user supplied bdsim options
    a.AddOptions(options)
    b.AddOptions(options)

    b.Write(outputfilename, overwrite=overwrite)
    if verbose:
        a.Write(outputfilename + "_raw", overwrite=overwrite)
        print('Total length: ',a.GetIntegratedLength())
        print('Total angle:  ',a.GetIntegratedAngle())
        print('items omitted: ')
        print(itemsomitted)
        print('number of omitted items: ',len(itemsomitted))

    return b,a,itemsomitted

def MadxTfs2GmadBeam(tfs, startname=None, verbose=False, gaussmatrix=False):
    print('Warning - using automatic generation of input beam distribution from madx tfs file - PLEASE CHECK!')

    if startname is None:
        startindex = 0
    else:
        try:
            startindex = tfs.IndexFromName(startname)
        except ValueError: # Then assume it's already an index.
            startindex = startname

    #MADX defines parameters at the end of elements so need to go 1 element
    #back if we can.

    if startindex > 0:
        startindex -= 1
    
    energy   = float(tfs.header['ENERGY'])
    gamma    = float(tfs.header['GAMMA'])
    particle = tfs.header['PARTICLE']
    ex       = tfs.header['EX']
    ey       = tfs.header['EY']
    sigmae   = float(tfs.header['SIGE'])
    sigmat   = float(tfs.header['SIGT'])

    if ex == 1:
        print('Horizontal emittance of 1 is too large - setting to 1e-9')
        ex = 1e-9
    if ey == 1:
        print('Horizontal emittance of 1 is too large - setting to 1e-9')
        ey = 1e-9
    
    data = tfs[startindex]

    if particle == 'ELECTRON' :
        particle = 'e-'
    elif particle == 'POSITRON' :
        particle = 'e+'
    elif particle == 'PROTON' :
        particle = 'proton'

    #print particle,energy,gamma,ex,ey
    if verbose:
        print('beta_x: ',data['BETX'],'alpha_x: ',data['ALFX'],'mu_x: ',data['MUX'])
        print('beta_y: ',data['BETY'],'alpha_y: ',data['ALFY'],'mu_y: ',data['MUY'])

    #gammax = (1.0+data['ALFX'])/data['BETX']
    #gammay = (1.0+data['ALFY'])/data['BETY']

    #note, in the main pybdsim.__init__.py Beam class is imported from Beam.py
    #so in this submodule when we do from .. import Beam it's actually the
    #already imported class that's being imported
    if gaussmatrix:
        beam   = _Beam.Beam(particle,energy,'gaussmatrix')
        for i, j in  [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                              (2, 2), (2, 3), (2, 4), (2, 5), (2, 6),
                                      (3, 3), (3, 4), (3, 5), (3, 6),
                                              (4, 4), (4, 5), (4, 6),
                                                      (5, 5), (5, 6),
                                                              (6, 6)]:
            beam["sigma{}{}".format(i, j)] = data["SIG{}{}".format(i, j)]
    else:
        beam   = _Beam.Beam(particle,energy,'gausstwiss')
        beam.SetBetaX(data['BETX'])
        beam.SetBetaY(data['BETY'])
        beam.SetAlphaX(data['ALFX'])
        beam.SetAlphaY(data['ALFY'])
        beam.SetDispX(data['DX'])
        beam.SetDispY(data['DY'])
        beam.SetDispXP(data['DPX'])
        beam.SetDispYP(data['DPY'])
        beam.SetEmittanceX(ex,'m')
        beam.SetEmittanceY(ey,'m')
        beam.SetSigmaE(sigmae)
    beam.SetXP0(data['PX'])
    beam.SetYP0(data['PY'])
    beam.SetX0(data['X'])
    beam.SetY0(data['Y'])
    
    return beam
