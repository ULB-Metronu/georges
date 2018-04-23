from re import split as _split

_allowedparticles = frozenset([
    'e-',
    'e+',
    'proton',
    'gamma'
])

_allowedflags = frozenset([
    '1',
    '2',
    '3'
])


class XSecBias(object):
    """
    A class for containing all information regarding cross section definitions.
    """
    def __init__(self, name, particle, processes, xsecfactors, flags):
        self.name        = self.SetName(name)
        self.particle    = self.SetParticle(particle)
        self.processlist = self.SetProcesses(processes)
        self.flaglist    = self.SetFlags(flags)
        self.xseclist    = self.SetXSecFactors(xsecfactors)
        self.CheckBiasedProcesses()

    def SetName(self, name):
        '''
        Set the bias name.
        Cannot be any upper/lowercase variant of reserved keyword "xsecBias".
        '''
        if name.lower() == 'xsecbias':
            raise ValueError("Forbidden name: "+ str(name) +".  Bias name cannot be any variant of 'xsecbias' (case insensitive)")
        return name

    def SetParticle(self, particle):
        """
        Set the particle for bias to be associated with.
        """
        if particle  in _allowedparticles:
            return particle
        else:
            raise ValueError("Unknown Particle type: " + str(particle) + ".")

    def SetProcesses(self, processes):
        '''
        Set the list of processes to be biased.
        processes hould be a space-delimited string of processes.
        '''
        processlist = _split(', |,| ', processes)
        if ('all' in processlist) and (len(processlist) > 1):
            raise ValueError("You cannot have both 'all' and another process be biased.")
        return processlist

    def SetFlags(self, flags):
        """
        Set flags.  flags should be a space-delimited string of integers, 1-3,
        in the same order as the processes,
        """
        flaglist = _split(', |,| ',flags)
        for number in flaglist:
            if number not in _allowedflags:
                raise ValueError("Unknown Flag type: " + str(number))
        return flaglist

    def SetXSecFactors(self, xsecs):
        """
        Set cross section factors.  xsecs should be a space-delimited string
        of floats, e.g. "1.0 1e13 1234.9"
        """
        xseclist = _split(', |,| ', xsecs)
        for factor in xseclist:
            if factor < 0:
                raise ValueError("Negative cross section factor: " +str(factor))
        return xseclist

    def CheckBiasedProcesses(self):
        if (len(self.xseclist) != len(self.flaglist)) or (len(self.xseclist) != len(self.processlist)):
            raise Warning("There must be a uniquely defined flag and xsecfactor to go  with every listed process.")

    def __repr__(self):
        particle = 'particle="' + self.particle + '"'
        proc     = 'proc="' + ' '.join(map(str, self.processlist)) + '"'
        xsec     = 'xsecfact={' + ','.join(map(str, self.xseclist)) + '}'
        flag     = 'flag={' + ','.join(map(str,self.flaglist)) + '}'

        s = self.name + ": xsecBias, " + ', '.join([particle,proc,xsec,flag]) + ';\n'

        return s
