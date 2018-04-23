from . import Data as _Data
import pylab as _pl

class Helper : 
    '''To help locate objects in the BDSIM visualisation, requires a BDSIM survey file'''
    
    def __init__(self, surveyFileName ) :
        self.survey      = _Data.Load(surveyFileName)
        self.x           = self.survey.X()
        self.y           = self.survey.Y()
        self.z           = self.survey.Z()
        self.worldCentre = self.getWorldCentre()
        

    def getWorldCentre(self, type = "linear") : 
        '''Returns the center in world coordinates of the centre of the visualisation space'''

        if type == "linear" : 
            xCentre = (self.x[0]+self.x[-1])/2.0
            yCentre = (self.y[0]+self.y[-1])/2.0
            zCentre = (self.z[0]+self.z[-1])/2.0
        elif type == "circular" : 
            xCentre = (self.x.max()+self.x.min())/2.0
            yCentre = (self.y.min()+self.y.max())/2.0
            zCentre = (self.z.max()+self.z.min())/2.0

       
        print("Visualisation.Helper.getWorldCentre>",xCentre,yCentre,zCentre)
        return _pl.array([xCentre,yCentre,zCentre])

    def findComponentCoords(self, componentName) :
        '''Returns the XYZ coordinates of a component relative to the centre'''
        names = self.survey.Name()

        matchingNames = []
        matchingIdx   = []

        idx = 0
        # loop over names 
        for name in names : 
            if name.rfind(componentName) != -1 : 
                print('Visualisation.Helper.findComponentCoords>',name,idx)
                matchingNames.append(name)
                matchingIdx.append(idx)
            idx = idx+1
            
        i = 0
        for match in matchingNames : 
            idx = matchingIdx[matchingNames.index(match)]
            x   = self.x[idx]
            y   = self.y[idx]
            z   = self.z[idx]
            v   = _pl.array([x,y,z]) 
            vp  = v - self.worldCentre
            print(match,idx,v,vp)
            
    def draw(self) : 
        '''Quick survey drawing for diagnostic reasons'''
        _pl.subplot(2,1,1);
        _pl.plot(self.z, self.x);
        _pl.plot(self.worldCentre[2], self.worldCentre[0],"+")
        _pl.xlabel("Z [m]")
        _pl.ylabel("X [m]")
        _pl.subplot(2,1,2);
        _pl.plot(self.z, self.y);
        _pl.plot(self.worldCentre[2], self.worldCentre[1],"+")
        _pl.xlabel("Z [m]")
        _pl.ylabel("Y [m]")


        
        pass
