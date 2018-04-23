import ROOT
from ROOT import TH1F
from ROOT import TFile

import glob
import optparse

def JoinRootHistograms(inputdir="./",outputfilename="output.root"):
    #get list of all root files in this dir
    files = glob.glob(inputdir+"*.root")

    print "Found ",len(files)," root files"

    hasextra = False #by default
    if len(files) == 0:
        print 'No root files found'
        return
    
    #open a single file as a template
    f = TFile(files[0])
    eloss = f.Get("ElossHisto")
    phits = f.Get("PhitsHisto")
    ploss = f.Get("PlossHisto")
    #test if it has extra histograms
    elosspe = f.Get("ElossPEHisto")
    if elosspe:
        hasextra = True
        phitspe = f.Get("PhitsPEHisto")
        plosspe = f.Get("PlossPEHisto")

    #detach the historams from the file so we don't loose
    #them when we close the file
    eloss.SetDirectory(0)
    phits.SetDirectory(0)
    ploss.SetDirectory(0)
    if hasextra:
        elosspe.SetDirectory(0)
        phitspe.SetDirectory(0)
        plosspe.SetDirectory(0)
        
    f.Close() #close template file

    #loop over rest of files and accumulate histograms
    if len(files) > 1:
        for filename in files[1:]:
            f = TFile(filename)
            elosstmp = f.Get("ElossHisto")
            phitstmp = f.Get("PhitsHisto")
            plosstmp = f.Get("PlossHisto")
            if hasextra:
                elosspetmp = f.Get("ElossPEHisto")
                phitspetmp = f.Get("PhitsPEHisto")
                plosspetmp = f.Get("PlossPEHisto")

            eloss += elosstmp
            phits += phitstmp
            ploss += plosstmp
            if hasextra:
                elosspe += elosspetmp
                phitspe += phitspetmp
                plosspe += plosspetmp
            f.Close()

    print "Analysed all files"
    print "Writing total histograms to ",outputfilename
    
    #open an output file for root
    f = TFile(outputfilename,"RECREATE")

    #write total histograms to root file
    f.WriteTObject(eloss)
    f.WriteTObject(phits)
    f.WriteTObject(ploss)
    if hasextra:
        f.WriteTObject(elosspe)
        f.WriteTObject(phitspe)
        f.WriteTObject(plosspe)

    f.Close()
    
#parser = optparse.OptionParser(usage)
#parser.add_option('-o','--outputfilename',action='store',default='output.root')
#parser.add_option('-d','--inputdir',action='store',default='./')

#options,args = parser.parse_args()

#if not options.outputfilename:
#    raise IOError("no output file specified")
#if not options.inputdir:
#    raise IOError("no input directory specified")
