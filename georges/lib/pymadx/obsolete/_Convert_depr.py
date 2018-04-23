import MadX as _MadX

def madx2gmad(iFilename, oFilename ) :
    o = open(oFilename,'w')
    
    t = _MadX.Twiss(iFilename)
    
    s = t.data['S'] 

    
    kw  = t.data['KEYWORD'] 
    na  = t.data['NAME'] 
    l   = t.data['L']
    k0l = t.data['K0L']
    k1l = t.data['K1L']
    k2l = t.data['K2L'] 

    i = 0 
    
    stw = str() 
    for v in s : 
        if kw[i] == 'MULTIPOLE' :
            if k0l[i] != 0 : 
                stw = na[i]+': dipole, l='+str(l[i])
            if k1l[i] != 0 : 
                stw = na[i]+': quadrupole, l='+str(l[i])+' k1='+str(k1l[i])
            if k2l[i] != 0 :
                stw = na[i]+': sextupole, l='+str(l[i])+' k1='+str(k1l[i])
        elif kw[i] == 'DRIFT' : 
            stw = na[i]+': drift, l='+str(l[i])

        stw += ';\n'
        
        o.write(stw)

        i += 1
