PDGind = {
#COMMON
    2212  : ('Proton',                'p'),
    -2212 : ('Antiproton',            r'$\={p}$'),
    2112  : ('Neutron',               'n'),
    -2112 : ('Neutron',               'n'),
    11    : ('Electron',              r'$e^{-}$'),
    -11   : ('Positron',              r'$e^{+}$'),
    22    : ('Photon',                r'$\gamma$'),
    12    : ('Electron Neutrino',     r'$\nu_{e}$'),
    -12   : ('Electron Antineutrino', r'$\=\nu_{e}$'),
    -13   : ('Antimuon',              r'$\mu^{+}$'),
#LIGHT MESONS I=1
    13    : ('Muon',                  r'$\mu^{-}$'),
    111   : ('Pion0',                 r'$\pi^{0}$'),
    211   : ('Pion+',                 r'$\pi^{+}$'),
    -211  : ('Pion-',                 r'$\pi^{-}$'),
#STRANGE MESONS
    321   : ('Kaon+',                 r'$K^{+}$'),
    -321  : ('Kaon-',                 r'$K^{-}$'),
    130   : ('K-Long',                r'$K_{L}^{0}$'),
    310   : ('K-Short',               r'$K_{S}^{0}$'),
#STRANGE BARYONS
    3122  : ('Lambda',                r'$\Lambda$'),
    3222  : ('Sigma+',                r'$\Sigma^{+}$'),
    3212  : ('Sigma0',                r'$\Sigma^{0}$'),
    3112  : ('Sigma-',                r'$\Sigma^{-}$'),
}

PDGname = {}
for k,v in PDGind.items():
    PDGname[v[0]] = k
    PDGname[v[0].lower()] = k
del k,v

def GetPDGInd(particlename):
    if particlename in PDGname:
        return PDGname[particlename]
    elif particlename.lower() in PDGname:
        return PDGname[particlename.lower()]
    else:
        raise ValueError("Unknown particle type")

def GetPDGName(particleid):
    try:
        return PDGind[particleid]
    except KeyError:
        print('Unknown particle id ',particleid)
        return ('','')
