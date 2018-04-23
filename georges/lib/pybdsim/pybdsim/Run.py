import os
import subprocess
import uuid
from . import Data

from . import _General


class ExecOptions(dict):
    def __init__(self,*args,**kwargs):
        """
        Executable options class for BDSIM. In addition, 'bdsimcommand' is an extra
        allowed keyword argument that allows the user to specify the bdsim executable
        of their choice.

        bdsimcommand='bdsim-devel'
        bdsimcommand='/Users/nevay/physics/bdsim-build2/bdsim'

        Based on python dictionary but with parameter checking.
        """
        self.__dict__.__init__()
        self._okFlags = ['batch',
                         'circular',
                         'generatePrimariesOnly',
                         'verbose']
        self._okArgs = ['bdsimcommand',
                        'file',
                        'output',
                        'outfile',
                        'ngenerate',
                        'seed',
                        'seedstate',
                        'survey',
                        'distrFile']
        for key,value in kwargs.iteritems():
            if key in self._okFlags or key in self._okArgs:
                self[key] = value
            else:
                raise ValueError(key+'='+str(value)+' is not a valid BDISM executable option')
        self.bdsimcommand = 'bdsim'
        if 'bdsimcommand' in self:
            self.bdsimcommand = self['bdsimcommand']
            self.pop("bdsimcommand", None)

    def GetExecFlags(self):
        result = dict((k,self[k]) for k in self.keys() if k in self._okFlags)
        return result

    def GetExecArgs(self):
        result = dict((k,self[k]) for k in self.keys() if k in self._okArgs)
        return result

class GmadModifier(object):
    def __init__(self, rootgmadfilename):
        self.rootgmadfilename = rootgmadfilename
        self.gmadfiles = [self.rootgmadfilename]
        self.DetermineIncludes(self.rootgmadfilename)
        self.CheckExtensions()

    def DetermineIncludes(self,filename) :
        f = open(filename)
        for l in f :
            if l.find("include") != -1 :
                includefile = l.split()[1].replace(";","")
                self.gmadfiles.append(includefile)
                self.DetermineIncludes(includefile)
        f.close()

    def CheckExtensions(self) :
        for filename in self.gmadfiles : 
            pass
    
    def ReplaceTokens(self,tokenDict) :
        pass
        


class Study(object):
    """
    A holder for multiple runs.
    """
    def __init__(self):
        self.execoptions = [] # exec options
        self.outputnames = [] # file names
        self.outputsizes = [] # in bytes

    def GetInfo(index=-1):
        """
        Get info about a particular run.
        """
        if index < 0:
            print("No runs yet")
            return
        
        i = index
        result = {'execoptions' : self.execoptions[i],
                  'outputname'  : self.outputnames[i],
                  'outputsize'  : self.outputsizes[i]}
        return result
                
    def Run(self, inputfile='optics.gmad',
            output='rootevent',
            outfile='output',
            ngenerate=1,
            bdsimcommand='bdsim-devel',
            **kwargs):
        eo = ExecOptions(file=inputfile, output=output, outfile=outfile, ngenerate=ngenerate, **kwargs)
        return self.RunExecOptions(eo)

    def RunExecOptions(self, execoptions, debug=False):
        if type(execoptions) != ExecOptions:
            raise ValueError("Not instance of ExecOptions")

        # shortcut
        eo = execoptions
        
        # prepare execution command
        command = eo.bdsimcommand
        for k in eo.GetExecFlags():
            command += ' --' + k
        for k,v in eo.GetExecArgs().iteritems():
            command += ' --' + str(k) + '=' + str(v)
        if debug:
            print('Command is')
            print(command)

        # send it to a log file
        outfilename = 'output'
        if 'outfile' in eo:
            outfilename = eo['outfile']
        command += ' > ' + outfilename + '.log'
        
        # execute process
        if debug:
            print('BDSIM Run')
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError:
            print('ERROR')
            return
        
        # get output file name - the latest file in the directory hopefully
        try:
            outfilename = eo['outfile']
        except KeyError:
            outfilename = _General.GetLatestFileFromDir(extension='*root') # this directory

        # record info
        self.execoptions.append(eo)
        self.outputnames.append(outfilename)
        try:
            self.outputsizes.append(os.path.getsize(outfilename))
        except OSError:
            self.outputsizes.append(0)


def RunBdsim(gmadpath, outfile, ngenerate=10000, batch=True,
             silent=False, options=None):
    """Runs bdsim with gmadpath as inputfile and outfile as outfile.
    Runs in batch mode by default, with 10,000 particles.  Any extra
    options should be provided as a string or iterable of strings of
    the form "--vis_debug" or "--vis_mac=vis.mac", etc.

    """
    args = ["bdsim",
            "--file={}".format(gmadpath),
            "--outfile={}".format(outfile),
            "--ngenerate={}".format(ngenerate)]
    if batch:
        args.append("--batch")

    if isinstance(options, basestring):
        args.append(options)
    elif options is not None:
        args.extend(options)

    if not silent:
        return subprocess.call(args)
    else:
        return subprocess.call(args, stdout=open(os.devnull, 'wb'))

def RunRebdsimOptics(rootpath, outpath, silent=False):
    """Run rebdsimOptics"""
    if not _General.IsROOTFile(rootpath):
        raise IOError("Not a ROOT file")
    if silent:
        return subprocess.call(["rebdsimOptics", rootpath, outpath],
                               stdout=open(os.devnull, 'wb'))
    else:
        return subprocess.call(["rebdsimOptics", rootpath, outpath])

def GetOpticsFromGMAD(gmad, keep_optics=False):
    """Get the optical functions as a BDSAsciiData instance from this
    GMAD file. If keep_optics is false then all intermediate files are
    discarded, otherwise the final optics ROOT file is written to ./
    """
    tmpdir = "/tmp/pybdsim-get-optics-{}/".format(uuid.uuid4())
    gmadname = os.path.splitext(os.path.basename(gmad))[0]
    os.mkdir(tmpdir)

    RunBdsim(gmad,
             "{}/{}".format(tmpdir, gmadname), silent=False,
             ngenerate=10000)
    bdsim_output_path = "{}/{}.root".format(tmpdir, gmadname)
    if keep_optics: # write output root file locally.
        RunRebdsimOptics(bdsim_output_path,
                         "./{}-optics.root".format(gmadname))
        return pybdsim.Data.Load("./{}-optics.root".format(gmadname))
    else: # do it in /tmp/
        RunRebdsimOptics(bdsim_output_path,
                         "{}/{}-optics.root".format(tmpdir, gmadname))
        return pybdsim.Data.Load("{}/{}-optics.root".format(tmpdir, gmadname))
