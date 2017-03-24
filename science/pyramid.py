#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
import numpy
import util.FITS
import base.aobase
class Pyramid(base.aobase.aobase):
    """
    pyramid wfs module.
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        if forGUISetup==1:
            npup=self.config.getVal("npup")
            print "TODO _ outputData size in pyramid"
            self.outputData=[(npup,npup),numpy.float32]
        else: # set up for simulation.
            self.npup=self.config.getVal("npup")
            obj=self.config.getVal("wfsOverview")
            wfsobj=obj.getWfsByID(idstr)
            self.pupil=wfsobj.pupil
            self.control={"parentType":self.parent.keys()[0]}
            self.nsteps=1#probably want more steps than 8 for larger modulation amplitudes
            if self.nsteps==0:
                self.nsteps=1
            pad=2
            self.pupAmp=numpy.zeros((self.npup*pad,self.npup*pad),numpy.complex64)
            self.tilt=numpy.zeros((self.nsteps,self.npup,self.npup),numpy.float32)
            tiltamp=0.5*self.npup/(pad*self.npup)*2*numpy.pi#this can be changed to change the modulation amplitude.
            
            xtiltfn=(numpy.fromfunction(lambda x,y:y,(self.npup,self.npup))-(self.npup-1)/2.)/self.npup
            tiltfn=-0.5*self.npup/float(self.npup*pad)*2*numpy.pi*(xtiltfn+xtiltfn.T)#this tilt centres it on 2x2 pixels.
            for i in range(self.nsteps):
                #first tilt everything so that its centred on 2x2 pixels (like we do with SHS images):
                self.tilt[i]=tiltfn
                #and now add some rotating tilt.
                if self.nsteps>1:
                    theta=i*2*numpy.pi/float(self.nsteps)+2*numpy.pi/(2*self.nsteps)
                    self.tilt[i]+=(tiltamp*xtiltfn*numpy.cos(theta)+tiltamp*xtiltfn.T*numpy.sin(theta))
            self.pyrImg=numpy.zeros((self.npup*pad,self.npup*pad),numpy.float32)
            #use pyrimg as temporary space to set up the pyramid phase tilt...
            s=self.npup*pad/2
            self.pyrImg[:s,:s]=0#(xtiltfn[:s,:s]*numpy.cos(numpy.pi/4)+xtiltfn[:s,:s].T*numpy.sin(numpy.pi/4))*0
            self.pyrImg[:s,s:]=xtiltfn.T*self.npup*numpy.pi#(xtiltfn[:s,:s]*numpy.cos(3*numpy.pi/4)+xtiltfn[:s,:s].T*numpy.sin(3*numpy.pi/4))*-10
            self.pyrImg[s:,s:]=(xtiltfn[:s,:s]+xtiltfn[:s,:s].T)*self.npup*numpy.pi
            self.pyrImg[s:,:s]=xtiltfn*self.npup*numpy.pi#(xtiltfn[:s,:s]*numpy.cos(5*numpy.pi/4)+xtiltfn[:s,:s].T*numpy.sin(5*numpy.pi/4))*0
            #now scale it
            self.pyrImg*=1.#change this to move the pupils in or out (larger moves them in)

            self.pyrPhaseMask=numpy.zeros((self.npup*pad,self.npup*pad),numpy.complex64)
            self.pyrPhaseMask.real=numpy.cos(self.pyrImg)
            self.pyrPhaseMask.imag=numpy.sin(self.pyrImg)

            self.outputData=self.pyrImg
            
    def calcPyrImg(self,phs):
        """Takes pupil plane, transports to focal plane (point of pyramid)
        scales for pixel scale, and then transports back to pupil plane.
        """
        self.pyrImg[:]=0
        for i in range(self.nsteps):
            self.pupAmp.real[:self.npup,:self.npup]=self.pupil*numpy.cos(phs+self.tilt[i])
            self.pupAmp.imag[:self.npup,:self.npup]=self.pupil*numpy.sin(phs+self.tilt[i])
            #go to the focus
            self.focAmp=(numpy.fft.fft2(self.pupAmp))
            #multiply by pyramid phase mask (shape of pyramid)
            self.focAmp*=self.pyrPhaseMask
            #and now transport back to pupil plane and detect
            pupilPlane=numpy.fft.ifft2(self.focAmp)
            self.pyrImg[:]+=numpy.abs(pupilPlane)**2

        return self.pyrImg

    def addNoise(self):
        return self.pyrImg#todo...
    
    def runPyramid(self,phs):
        self.calcPyrImg(phs)
        self.addNoise()
        
    def generateNext(self,ms=None):
        current=self.control["parentType"]
        if self.generate==1:
            if self.parent[current].dataValid==0:
                print("INFORMATION:pyramid: waiting for data from dm")
                self.dataValid=0
            else:
                self.runPyramid(self.parent[current].outputData)
                self.dataValid=1
        else:
            self.dataValid=0

        
def test():
    import util.tel
    class Dummy:
        dataValid=1
        outputData=numpy.zeros((160,160),"f")
        pupil=util.tel.Pupil(160,80,0).fn
        
    class Overview:
        def getWfsByID(self,name):
            return Dummy()

    class Config:
        def __init__(self):
            self.dict={"npup":160,
                       "wfsOverview":Overview(),
            }
            self.rank=0
        def setSearchOrder(self,o=None):
            pass
        def getVal(self,name,default=None,raiseerror=0):
            if self.dict.has_key(name):
                return self.dict[name]
            elif raiseerror==0:
                return default
            elif default is None:
                raise Exception("%s not found"%name)
            else:
                return default
        
    parent=Dummy()
    config=Config()
    p=Pyramid({1:parent},config)
    return p

if __name__=="__main__":
    p=test()
    p.generateNext()
