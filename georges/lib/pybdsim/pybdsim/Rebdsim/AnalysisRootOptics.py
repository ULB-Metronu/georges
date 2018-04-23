import AnalysisRoot 

class AnalysisRootOptics : 
    def __init__(self, filelist ) : 
        self.ar = AnalysisRoot.AnalysisRoot(filelist) 
        
    def LoopOverSamplers(self) :         
        for s in self.ar.samplerNames : 
            sd = self.SamplerAnalysis(self.ar.samplerDict[s])
            print s,sd

    def SamplerAnalysis(self, chain) :         
        chain.SetBranchStatus("x",1) 
        chain.SetBranchStatus("y",1) 
        chain.SetBranchStatus("z",1) 
        
        xmean = 0.0 
        ymean = 0.0 
        zmean = 0.0 
        nevt  = 0

        for e in range(0,chain.GetEntries(),1) : 
            chain.GetEntry(e)
            xmean += chain.x 
            ymean += chain.y
            zmean += chain.z
            nevt  += 1 
#            print chain.x, chain.y, chain.z 
        xmean /= nevt
        ymean /= nevt
        zmean /= nevt

        return [xmean,ymean,zmean]



    
    
