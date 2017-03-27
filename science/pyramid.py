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
        if type(parent)!=type({}):
            parent={"closed":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.npup=self.config.getVal("npup")
        obj=self.config.getVal("wfsOverview")
        wfsobj=obj.getWfsByID(idstr)
        self.nimg=wfsobj.nimg
        self.imageOnly=config.getVal("imageOnly",default=0)
        
        if forGUISetup==1:
            if self.imageOnly==0:
                self.outputData=[(self.nimg/2,self.nimg/2,2),numpy.float32]
            else:
                self.outputData=[(self.nimg,self.nimg),numpy.float32]
        else: # set up for simulation.
            self.pupil=wfsobj.pupil
            self.nfft=wfsobj.nfft
            self.clipsize=wfsobj.clipsize
            self.nimg=wfsobj.nimg
            if self.imageOnly==0:
                self.outputData=numpy.zeros((self.nimg/2,self.nimg/2,2),numpy.float32)
            else:
                self.outputData=numpy.zeros((self.nimg,self.nimg),numpy.float32)
            self.control={"parentType":self.parent.keys()[0]}
            self.nsteps=wfsobj.pyrSteps#probably want more steps than 8 for larger modulation amplitudes
            if self.nsteps==0:
                self.nsteps=1
            self.pupAmp=numpy.zeros((self.nfft,self.nfft),numpy.complex64)
            self.tilt=numpy.zeros((self.nsteps,self.npup,self.npup),numpy.float32)
            tiltamp=0.5*self.npup/(self.nfft)*2*numpy.pi*wfsobj.pyrModAmp#this can be changed to change the modulation amplitude.
            
            xtiltfn=(numpy.fromfunction(lambda x,y:y,(self.npup,self.npup))-(self.npup-1)/2.)/self.npup
            tiltfn=-0.5*self.npup/float(self.nfft)*2*numpy.pi*(xtiltfn+xtiltfn.T)#this tilt centres it on 2x2 pixels.
            for i in range(self.nsteps):
                #first tilt everything so that its centred on 2x2 pixels:
                self.tilt[i]=tiltfn# (like we do with SHS images)
                #and now add some rotating tilt.
                if self.nsteps>1:
                    theta=i*2*numpy.pi/float(self.nsteps)+2*numpy.pi/(2*self.nsteps)
                    self.tilt[i]+=(tiltamp*xtiltfn*numpy.cos(theta)+tiltamp*xtiltfn.T*numpy.sin(theta))
            self.pyrImg=numpy.zeros((self.nfft,self.nfft),numpy.float32)
            #use pyrimg as temporary space to set up the pyramid phase tilt...
            s=self.nfft/2
            xtiltfn=(numpy.fromfunction(lambda x,y:y,(s,s))-(s-1)/2.)/s
            self.pyrImg[:s,:s]=(xtiltfn[:s,:s]+xtiltfn[:s,:s].T)
            self.pyrImg[:s,s:]=(-xtiltfn[:s,:s]+xtiltfn[:s,:s].T)
            self.pyrImg[s:,s:]=(-xtiltfn[:s,:s]-xtiltfn[:s,:s].T)
            self.pyrImg[s:,:s]=(xtiltfn[:s,:s]-xtiltfn[:s,:s].T)
            #now scale it: change this to move the pupils in or out
            #By 0.5 moves in by half a diameter, 1.5 out by half a diameter.
            self.pyrImg*=numpy.pi/2*self.npup#larger=move out

            
            self.pyrPhaseMask=numpy.zeros((self.nfft,self.nfft),numpy.complex64)
            self.pyrPhaseMask.real=numpy.cos(self.pyrImg)
            self.pyrPhaseMask.imag=numpy.sin(self.pyrImg)

            #self.outputData=self.pyrImg
            self.takeRefSlopes()

            self.sentPlotsCnt=0
            
    def calcPyrImg(self,phs):
        """Takes pupil plane, transports to focal plane (point of pyramid)
        scales for pixel scale, and then transports back to pupil plane.
        """
        self.pyrImg[:]=0
        s=(self.nfft-self.npup)//2
        for i in range(self.nsteps):
            self.pupAmp.real[s:s+self.npup,s:s+self.npup]=self.pupil*numpy.cos(phs+self.tilt[i])
            self.pupAmp.imag[s:s+self.npup,s:s+self.npup]=self.pupil*numpy.sin(phs+self.tilt[i])
            #go to the focus
            self.focAmp=(numpy.fft.fft2(self.pupAmp))
            #multiply by pyramid phase mask (shape of pyramid)
            self.focAmp*=self.pyrPhaseMask
            #and now transport back to pupil plane and detect
            pupilPlane=numpy.fft.ifft2(self.focAmp)
            self.pyrImg[:]+=numpy.abs(pupilPlane)**2

        return self.pyrImg

    def binPyr(self):#bin down onto the detector
        s=(self.nfft-self.clipsize)/2
        img=self.pyrImg[s:s+self.clipsize,s:s+self.clipsize]
        binfact=self.clipsize/self.nimg
        if self.clipsize%self.nimg!=0:
            raise Exception("Can only do integer binning at present")
        img.shape=self.nimg,binfact,self.nimg,binfact
        img=img.sum(3).sum(1)
        self.binImg=img
        return self.binImg
    
    def addNoise(self):#readout the detector
        print "TODO: Add pyramid noise"
        return self.binImg#todo...

    def takeRefSlopes(self):
        if self.imageOnly==0:
            self.refSlopes=None
            self.calcPyrImg(0)
            self.binPyr()
            self.computeSlopes()
            self.refSlopes=self.outputData.copy()
    
    def computeSlopes(self):
        s=self.binImg.shape[0]/2
        img=self.binImg
        tot=(img[:s,:s]+img[s:,:s]+img[:s,s:]+img[s:,s:])/4.
        x=(img[:s,:s]+img[s:,:s]-img[:s,s:]-img[s:,s:])/tot
        y=(img[:s,:s]-img[s:,:s]+img[:s,s:]-img[s:,s:])/tot
        if self.refSlopes is not None:
            x-=self.refSlopes[:,:,0]
            y-=self.refSlopes[:,:,1]
        self.outputData[:,:,0]=x
        self.outputData[:,:,1]=y
        return x,y
    
    def runPyramid(self,phs):
        self.calcPyrImg(phs)
        self.binPyr()
        self.addNoise()
        if self.imageOnly==0:
            self.computeSlopes()
        else:
            self.outputData[:]=self.binImg
            
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

    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        txt=""
        if self.idstr==None or self.idstr=="":
            id=""
        else:
            id=" (%s)"%self.idstr
        if self.sentPlotsCnt==0:
            #outputData is only valid for one object at a time, when that has just run...
            txt+="""<plot title="Pyr img (noiseless) %s" cmd="data=%s.pyrImg" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            txt+="""<plot title="Pyr binned img %s" cmd="data=%s.binImg" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            txt+="""<plot title="X slopes %s" cmd="data=%s.outputData[:,:,0]" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            txt+="""<plot title="Y slopes %s" cmd="data=%s.outputData[:,:,1]" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
        return txt
        
def test():
    import util.tel
    import util.dm
    class Dummy:
        dataValid=1
        outputData=numpy.zeros((160,160),"f")
        #outputData[:]=numpy.arange(160)/160.
        ms=util.dm.MirrorSurface("pspline",160,11)
        actmap=numpy.zeros((11,11),"f")
        actmap[3,3]=1
        ms.fit(actmap)
        outputData[:]=ms.phsOut
        pupil=util.tel.Pupil(160,80,0).fn
        npup=160
        nfft=npup*2
        clipsize=nfft
        nimg=20#clipsize/2#(nact-1)*2
        
    class Overview:
        def getWfsByID(self,name):
            return Dummy()

    class Config:
        def __init__(self):
            self.dict={"npup":160,
                       "wfsOverview":Overview(),
                       "imageOnly":0,
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
    centx=p.x
    centy=p.y
