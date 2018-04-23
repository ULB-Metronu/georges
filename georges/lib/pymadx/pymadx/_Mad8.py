# http://mad.web.cern.ch/mad/madx.old/mad8tomadX/
# constant 
# multipole

import os as _os 
import re as _re 

def ConvertDir(inputDir) : 
    files = _os.listdir(inputDir) 
    for f in files : 
        if f.find('~') == -1 : 
            Mad8ToMadX(inputDir+'/'+f)
        
def Mad8ToMadX(inputName) : 
    inputFile  = open(inputName) 
    outputName  = inputName[inputName.rfind('/')+1:inputName.rfind('.')]+".xsifx"
    outputFile = open(outputName,'w')
    
    print 'Mad8ToMadX, input > ',inputName
    print 'Mad8ToMadX, output> ',outputName


    ml = '' # merged line 
    inputFile1 = []
    for l in inputFile : 
        l  = l.rstrip() 
        
        # find line continuations (apersand index) 
        ai = l.find('&')
        if ai != -1 : 
            l = l.replace('&','')
            ml = ml+l
        else :
            if len(ml) == 0 :
                ml = l
            else :
                ml = ml+l
            inputFile1.append(ml)
            ml = ''
        
    ltc = ';'
    for l in inputFile1 : 
        l = l.rstrip()
        
        # find comments (comment index)
        ci = l.find('!')

        # Deal with each line
        if ci == -1 and len(l) == 0: # Empty line  
            pl = l
        elif ci == -1 and len(l) > 0 : # No comment and actual line 
            pl = l+ltc
        elif ci == 0 : # Whole line is comment
            pl = l 
        else : # comment present 
            if len(l) == 0 : 
                print ci
            # split line on comment 
            bcl = l[0:ci] # before comment
            acl = l[ci:]  # after comment
            pl  = bcl+ltc+acl
 
        #######################################################################################
        # Simple replacements 
        #######################################################################################
        # replace intr with intrument 
        if pl.find('INSTRUMENT') == -1 :
            pl = pl.replace('INST','INSTRUMENT')        
        # replace moni with monitor
        if pl.find('MONITOR') == -1 :
            pl = pl.replace('MONI','MONITOR')
        # replace mark with marker 
        if pl.find('MARKER') == -1 :
            pl = pl.replace('MARK','MARKER')
        # replace drif with drift 
        if pl.find('DRIFT') == -1 :
            pl = pl.replace('DRIF','DRIFT')
        # replace quad with quadrupole 
        if pl.find('QUADRUPOLE') == -1 :
            pl = pl.replace('QUAD','QUADRUPOLE') 
        # replace sext with sextupole 
        if pl.find('SEXTUPOLE') == -1 : 
            pl = pl.replace('SEXT','SEXTUPOLE')
        # replace octu with octupole
        if pl.find('OCTUPOLE') == -1 : 
            pl = pl.replace('OCTU','OCTUPOLE')
        # replace sben with sbend
        if pl.find('SBEND') == -1 :
            pl = pl.replace('SBEN','SBEND')
        # replace aper with aperture 
        if pl.find('APERTURE') == -1 :
            pl = pl.replace('APER','APERTURE') 
        # replce ecoll with ecollimator
        if pl.find('ECOLLIMATOR') == -1 :
            pl = pl.replace('ECOLL','ECOLLIMATOR') 
        if pl.find('ECOLLIMATOR') == -1 :
            pl = pl.replace('ECOL','ECOLLIMATOR') 
        # replce rcoll with ecollimator
        if pl.find('RCOLLIMATOR') == -1 :
            pl = pl.replace('RCOLL','RCOLLIMATOR') 
        if pl.find('RCOLLIMATOR') == -1 :
            pl = pl.replace('RCOL','RCOLLIMATOR')
        if pl.find('HKICKER') == -1 : 
            pl = pl.replace('HKICK','HKICKER')
        if pl.find('HKICKER') == -1 : 
            pl = pl.replace('HKIC','HKICKER')
        if pl.find('VKICKER') == -1 : 
            pl = pl.replace('VKICK','HKICKER')
        if pl.find('VKICKER') == -1 : 
            pl = pl.replace('VKIC','HKICKER')



        #######################################################################################
        # Regular expressions 
        #######################################################################################
        # constant 
        if pl.find('CONSTANT') != -1 : 
#            print pl
            m = _re.search('(\w+)\s*:\s*CONSTANT\s*=\s*([A-Za-z0-9+-/*.()]+)',pl)
            pl = '   const '+m.group(1)+'='+m.group(2)+';'
            print pl

        # multipole 
        if pl.find('MULTIPOLE') != -1: 
            print pl
            n  = _re.search('([A-Za-z0-9.]+)\s*:\s*MULTIPOLE',pl)
            m0 = _re.search('K0L\s*=\s*([A-Za-z0-9+-/*.()]+),',pl)
            m1 = _re.search('K1L\s*=\s*([A-Za-z0-9+-/*.()]+),',pl)
            m2 = _re.search('K2L\s*=\s*([A-Za-z0-9+-/*.()]+),',pl)
            m3 = _re.search('K3L\s*=\s*([A-Za-z0-9+-/*.()]+),',pl)
            m4 = _re.search('K4L\s*=\s*([A-Za-z0-9+-/*.()]+),',pl)
            m5 = _re.search('K5L\s*=\s*([A-Za-z0-9+-/*.()]+),',pl)
            m6 = _re.search('K6L\s*=\s*([A-Za-z0-9+-/*.()]+),',pl)
            mt0 = _re.search('T0\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)       # match t0
            mt1 = _re.search('T1\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)       # match t1
            mt2 = _re.search('T2\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)       # match t2
            mt3 = _re.search('T3\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)       # match t3
            mt4 = _re.search('T4\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)       # match t4
            mt5 = _re.search('T5\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)       # match t5
            mt6 = _re.search('T6\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)       # match t6
            ml  = _re.search(' L\s*=\s*([0-9.E+-])',pl)                    # match length

            mt  = _re.search('TYPE\s*=\s*("(.*?)")',pl)                    # match type
            mr  = _re.search('LRAD\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl)     # match lrad
            ma  = _re.search('APERTURE\s*=\s*([A-Za-z0-9+-/*.()]+)\s*',pl) # match aperture

            pl = n.group(1)+' : MULTIPOLE '

            if ml : 
                pl = pl+', L='+ml.group(1)

            k0l = '0.0'; t0 = '0.0'
            k1l = '0.0'; t1 = '0.0'
            k2l = '0.0'; t2 = '0.0'
            k3l = '0.0'; t3 = '0.0'
            k4l = '0.0'; t4 = '0.0'
            k5l = '0.0'; t5 = '0.0'
            k6l = '0.0'; t6 = '0.0'
                
            if m0 : 
                k0l = m0.group(1)
            if m1 : 
                k1l = m1.group(1)
            if m2 : 
                k2l = m2.group(1)
            if m3 : 
                k3l = m3.group(1)
            if m4 : 
                k4l = m4.group(1)
            if m5 : 
                k5l = m5.group(1)
            if m6 : 
                k6l = m6.group(1)

            if mt0 : 
                t0 = mt0.group(1)
            if mt1 : 
                t1 = mt1.group(1)
            if mt2 : 
                t2 = mt2.group(1)
            if mt3 : 
                t3 = mt3.group(1)
            if mt4 : 
                t4 = mt4.group(1)
            if mt5 : 
                t5 = mt5.group(1)
            if mt6 : 
                t6 = mt6.group(1)            

            pl = pl+ ', KNL={'+k0l+','+k1l+','+k2l+','+k3l+','+k4l+','+k5l+','+k6l+'}'
            
            if mt : 
                pl = pl +', TYPE='+mt.group(1)
            if ma : 
                pl = pl +', APERTURE='+ma.group(1)
            if mr : 
                pl = pl +', LRAD='+mr.group(1)

            pl = pl+';'

            print pl
                
        # attributes 
        def maFunc(m) : 
            return m.group(1)+'->'+m.group(2)

        pl = _re.sub('([a-zA-Z0-9_-]+)\[([a-zA-Z]+)\]',maFunc,pl) # replace [] with ->
        outputFile.write(pl+'\n')
