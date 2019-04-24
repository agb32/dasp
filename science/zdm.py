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
#This used to be glao_zdm, but name was changed as it works for anything...
import numpy
import util.zernikeMod
import util.FITS
#from cmod.interp import gslCubSplineInterp
import base.aobase
from scipy.signal import lfilter

class dm(base.aobase.aobase):
    """
    Zernike DM simulation object:
    Perfect Zernike fitter.
    Add DM phases to input phase array and output
    Update mirror figure on request from reconstructor
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        if type(parent)!=type({}):
            parent={"1":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.sendFullDM=self.config.getVal("sendFullDM",default=0)#used if connecting to wideField.py science module
        self.atmosPhaseType=self.config.getVal("atmosPhaseType",default="phaseonly")
        if forGUISetup==1:
            self.dmObj=self.config.getVal("dmOverview",default=self.config.getVal("dmObj"),raiseerror=0)
            if self.dmObj.getDM(idstr).sendFullDM:
                dmpup=self.dmObj.calcdmpup(self.idstr[0])#number of pixels to store the phase. May be >npup if not ground conjugate.
                if self.atmosPhaseType=="phaseonly":
                    self.outputData=[(dmpup,dmpup),numpy.float32]
                else:
                    self.outputData=[(2,dmpup,dmpup),numpy.float32]
            else:
                npup=self.config.getVal("npup")
                if self.atmosPhaseType=="phaseonly":
                    self.outputData=[(npup,npup),numpy.float32]
                else:
                    self.outputData=[(2,npup,npup),numpy.float32]
        else: # set up for simulation.
            self.npup=self.config.getVal("npup")
            self.atmosGeom=self.config.getVal("atmosGeom",default=None,raiseerror=0)
            self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
            if self.dmObj==None or type(self.dmObj)!=type(self.atmosGeom):
                print "DEPRECATION: warning: dmObj should now be dmOverview"
                self.dmObj=self.config.getVal("dmObj",default=None,raiseerror=0)
            self.interpolationNthreads=self.config.getVal("interpolationNthreads",default=0)
            self.dmInfo=self.dmObj.getDM(self.idstr[0])
            self.monteNoll=None
            self.montePhaseCovfname=self.config.getVal("montePhaseCovFilename",raiseerror=0)
            self.actStart=None#offsets into the parent actuator lists...
            self.actEnd=None
            self.sourceLam=self.dmInfo.reconLam
	    self.dmpup=self.dmObj.calcdmpup(self.idstr[0])#number of pixels to store the phase. May be >npup if not ground conjugate.
            self.pupil=self.config.getVal("pupil",raiseerror=0)
            if self.pupil==None:
                r2=0
            else:
                r2=self.pupil.r2
            self.dmpupil=self.dmInfo.computeDMPupil(self.atmosGeom,r2)[2]
	    self.conjHeight=self.dmObj.getHeight(self.idstr[0])
	    self.actoffset=self.dmObj.getactoffset(self.idstr[0])

            self.nmodes=self.dmInfo.nact#number of zernikes to use.
            if self.dmInfo.zonalDM==1:
                raise Exception("zdm does not work for a zonal DM (label %s) - must be modal"%self.dmInfo.label)
            self.dmphs=numpy.zeros((self.dmpup,self.dmpup),numpy.float32)            # DM figure
            self.subpxlInterp=self.dmInfo.subpxlInterp
            if self.conjHeight==0:
                self.subpxlInterp=0
            if self.subpxlInterp:
                #self.interpolated=numpy.zeros((self.npup,self.npup),numpy.float32)
                self.yaxisInterp=numpy.arange(self.npup+1).astype(numpy.float64)
                self.xaxisInterp=numpy.arange(self.npup+1).astype(numpy.float64)
            else:
                #self.interpolated=None
                self.yaxisInterp=None
                self.xaxisInterp=None
            #self.gamma=config.getVal("gamma")#gain for zernikes.  Can be a float or an array length nmodes.
            #There are 2 ways of getting zernikes - the original, using RWW way, which returns circular zernikes, or FA way, which allows a pupil mask to be specified (and so zernikes aren't always orthogonal)
            #self.sourceID=self.dmObj.getSourceID(self.idstr[0])
            #self.sourceLam=self.atmosGeom.sourceLambda(self.sourceID)#self.config.getVal("sourceLam")#wavelength for which the mirror will be shaped.
            #self.wfslam=self.config.getVal("wfslam")#reconstructor wavelength
	    self.zoffset=self.config.getVal("zoffset",default=None,raiseerror=0)#zernike offsets that can be applied.  An array of length nmodes, ie numpy.zeros((self.nmodes,),numpy.float64)
            self.orthonormalZernike=self.config.getVal("orthonormalZernike",default=1)#DONT set this for MAP.

	    self.reconData=None
            #self.wavelengthRatio=self.wfslam/self.sourcelam
            self.telDiam=self.config.getVal("telDiam")
            if self.sendFullDM:
                self.outputData=self.dmphs
                if self.atmosPhaseType!="phaseonly":
                    print "Warning - sendFullDM selected for phase type %s.  May not work..."%self.atmosPhaseType
            else:
                if self.atmosPhaseType=="phaseonly":
                    self.outputData=numpy.zeros((self.npup,self.npup),numpy.float32)
                else:
                    self.outputData=numpy.zeros((2,self.npup,self.npup),numpy.float32)
                    self.outputData[1]=1.#set to 1... in case its a poking simulation.

            # Make Zernike fns over DM, and
            slow=0
            if slow:#not as accurate either I think.
                from cmod.zernike import zern
                import math
                self.zern=numpy.zeros((self.nmodes,self.dmpup,self.dmpup),numpy.float64)             # make the Zernikes over the pupil
                self.nzrad=util.zernikeMod.nm(self.nmodes)[0]#config.getVal("nzrad")

                pp=[0]
                qq=[0]
                tt=[1]
                trig=0
                for p in range(1,self.nzrad+1):
                    for q in range(p+1):
                        if(math.fmod(p-q,2)==0):
                            if(q>0):
                                pp.append(p)
                                qq.append(q)
                                trig=1-trig
                                tt.append(trig)
                                pp.append(p)
                                qq.append(q)
                                trig=1-trig
                                tt.append(trig)
                            else:
                                pp.append(p)
                                qq.append(q)
                                tt.append(1)
                                trig=1-trig
                for j in range(self.nmodes):
                    ztmp=numpy.zeros((self.dmpup,self.dmpup),numpy.float64)
                    zern(ztmp,pp[j],qq[j],tt[j])
                    self.zern[j]=ztmp
                #util.FITS.Write(self.zern,'zernike_dump.fits')#agb: bother doing this?
            else:#use the faster method.
                self.zern=util.zernikeMod.Zernike(self.dmpupil,self.nmodes,computeInv=0).zern
                if self.orthonormalZernike==1:#DON'T SET THIS FOR MAP
                    print "Scaling zernikes (don't do this for MAP)"
                    util.zernikeMod.normalise(self.zern)#normalise to orthonormal (ie numpy.sum(zern[i]*zern[i])==1).
                elif self.orthonormalZernike==2:#scale for NOLL:
                    print "Scaling zernikes for NOLL"
                    util.zernikeMod.normalise(self.zern,scaleto=self.zern[0].sum()**2)
            self.control={"dm_update":1,"phaseCovariance":0}
            self.lastPhaseCovariance=0

            for i in xrange(len(self.idstr)):
                self.initialise(self.parentList[i],self.idstr[i])


    def initialise(self,parentDict,idstr):
        """note parent may be a dictionary here (atmos and recon)"""
        if type(parentDict)!=type({}):
            parentDict={"1":parentDict}
        this=base.aobase.resourceSharer(parentDict,self.config,idstr,self.moduleName)
        npup=self.npup
        self.thisObjList.append(this)

        this.sourceID=self.dmObj.getSourceID(idstr)
        this.sourceAlt=self.atmosGeom.sourceAlt(this.sourceID)
        this.sourceLam=self.atmosGeom.sourceLambda(this.sourceID)
        if this.sourceLam==None:
            this.sourceLam=this.config.getVal("sourceLam")
            print "WARNING - DM sourceLam not found in atmosGeom, using %d"%this.sourceLam
        this.sourceTheta=self.atmosGeom.sourceTheta(this.sourceID)*numpy.pi/180/3600.
        this.sourcePhi=self.atmosGeom.sourcePhi(this.sourceID)*numpy.pi/180

        wavelengthAdjustor=self.sourceLam/this.sourceLam#this.wavelengthRat
        
        this.lineOfSight=util.dm.DMLineOfSight(self.dmpup,self.npup,self.conjHeight,self.dmphs,this.sourceAlt,this.sourceTheta,this.sourcePhi,self.telDiam,wavelengthAdjustor,self.xaxisInterp,self.yaxisInterp,dmTiltAngle=0.,dmTiltTheta=0.,alignmentOffset=(0,0),subpxlInterp=self.subpxlInterp,pupil=self.pupil,nthreads=self.interpolationNthreads)
        """
        #for a source with theta and a conjugate height c, the
        #separation at this height will be c tan(theta).  As a
        #fraction of the phasescreen width, this will be c
        #tan(theta)/telDiam.  So, in pixels, this will be c
        #tan(theta)/telDiam * npup.  Multiply this by
        #cos/sin(phi), and you then have the correct offsets.
        this.xoff=(self.dmpup-self.npup)/2#central in the case of
        this.yoff=(self.dmpup-self.npup)/2#ground conjugated.
        this.xoffend=this.xoff+npup
        this.yoffend=this.xoff+npup

        this.sourceID=self.dmObj.getSourceID(idstr)
        this.sourceAlt=self.atmosGeom.sourceAlt(this.sourceID)

        if self.conjHeight!=0:
            this.sourceTheta=self.atmosGeom.sourceTheta(this.sourceID)*numpy.pi/180/3600.
            this.sourcePhi=self.atmosGeom.sourcePhi(this.sourceID)*numpy.pi/180
            if this.sourceAlt>0 and self.subpxlInterp==0:
                raise Exception("zdm - finite altitude source specified but subpxlInterp==0.")
            
            r=self.conjHeight*numpy.tan(this.sourceTheta)/self.telDiam*self.npup
            if self.subpxlInterp:
                if this.sourceAlt>0:
                    if this.sourceAlt<self.conjHeight:
                        print "zdm Warning - you have a DM conjugated above a source height - is this what you expect (actuators will be reversed)"
                    x=r*numpy.cos(this.sourcePhi)
                    y=r*numpy.sin(this.sourcePhi)
                    #xoff+x would be the starting point if source at infinity.
                    #xoff+x+npup would be the ending point.
                    #So, I need to reduce this by a factor of (alt-conj)/alt
                    npxl=self.npup*(this.sourceAlt-self.conjHeight)/this.sourceAlt
                    #if npxl <0, means have to flip the mirror...
                    this.xoff=(this.xoff+x)+(npup-numpy.fabs(npxl))/2.
                    this.yoff=(this.yoff+y)+(npup-numpy.fabs(npxl))/2.
                    this.xoffsub=this.xoff%1
                    this.yoffsub=this.yoff%1
                    tol=1e-10#allow for floating point inprecision 
                    this.xoffend=int(numpy.ceil(this.xoff+numpy.fabs(npxl)-tol))#oversize by 1pxl
                    this.yoffend=int(numpy.ceil(this.yoff+numpy.fabs(npxl)-tol))#unless exact fit
                    this.xoff=int(this.xoff+tol)
                    this.yoff=int(this.yoff+tol)
                    if npxl>0:#source above DM conjugate...
                        this.xaxisInterp=numpy.arange(self.npup).astype(numpy.float64)*(npxl-1.)/(self.npup-1.)+this.xoffsub
                        this.yaxisInterp=numpy.arange(self.npup).astype(numpy.float64)*(npxl-1.)/(self.npup-1.)+this.yoffsub
                    else:#need to flip since will see mirror in reverse (I think)
                        this.xaxisInterp=numpy.arange(self.npup-1,-1,-1).astype(numpy.float64)*(numpy.fabs(npxl)-1.)/(self.npup-1.)+this.xoffsub
                        this.yaxisInterp=numpy.arange(self.npup-1,-1,-1).astype(numpy.float64)*(numpy.fabs(npxl)-1.)/(self.npup-1.)+this.yoffsub
                else:#source at infinity.
                    this.xoff+=r*numpy.cos(this.sourcePhi)
                    this.yoff+=r*numpy.sin(this.sourcePhi)
                    this.xoffsub=this.xoff%1
                    this.yoffsub=this.yoff%1
                    this.xoff=int(this.xoff)
                    this.yoff=int(this.yoff)
                    ax=ay=1
                    if this.yoffsub==0:
                        ay=0
                    if this.xoffsub==0:
                        ax=0
                    this.yoffend=this.yoff+self.npup+ay#oversize by 1 pupil unless exact fit.
                    this.xoffend=this.xoff+self.npup+ax
                    this.xaxisInterp=numpy.arange(self.npup).astype(numpy.float64)+this.xoffsub
                    this.yaxisInterp=numpy.arange(self.npup).astype(numpy.float64)+this.yoffsub
            else:
                this.xoff+=int(r*numpy.cos(this.sourcePhi)+0.5)#round
                this.yoff+=int(r*numpy.sin(this.sourcePhi)+0.5)#correctly
                this.xoffend=this.xoff+self.npup
                this.yoffend=this.yoff+self.npup
                    

            
        # select the correct part of the dm phs.
        if this.xoff<0 or this.yoff<0 or this.xoffend>self.dmpup or this.yoffend>self.dmpup:
            raise Exception("DM pupil not large enough to hold all the conjugated sources... %d %d %d %d %d"%(this.xoff,this.yoff,self.xoffend,self.yoffend,self.dmpup))

        this.selectedDmPhs=self.dmphs[this.yoff:this.yoffend,this.xoff:this.xoffend]

        this.sourceLam=self.atmosGeom.sourceLambda(this.sourceID)#this.config.getVal("sourceLam")          # Wavelength for this optical path
        #this.wfsLam=this.config.getVal("wfslam",default=this.sourceLam)
        #this.wavelengthRatio=this.wfsLam/this.sourceLam#wfs_lam/lam
        this.wavelengthAdjustor=self.sourceLam/this.sourceLam#this.wavelengthRatio/self.wavelengthRatio#the mirror will be shaped as for self.wavelengthRatio... if this.sourceLam is longer than self.sourceLam, radians of phase P-V will be smaller, so less change needed in wavelengthAdjustor.
        """

    def finalInitialisation(self):
        """Just check that there aren't 2 objects with same idstr..."""
        tmp=[]
        for id in self.idstr:
            if id in tmp:
                raise Exception("xinterp_dm - resource sharing cannot have 2 objects with same idstr")
            tmp.append(id)

    def addNewIdObject(self,parent,idstr):
        """Add a new ID string - for resource sharing
        This should be overridden for objects that allow resource sharing.  For objects that don't, it shouldn't be used.
        Overriding this method should include all that is necessary to read the config file, and get all necessary info set up for this particular idstr.
        """
        self.idstr.append(idstr)
        if type(parent)!=type({}):
            parent={"1":parent}
        self.parentList.append(parent)
        self.initialise(parent,idstr)
        return self

    def newParent(self,parent,idstr=None):
        """Most likely - if created by simsetup GUI, the object will
        originally have None as parent, and newParent will be called
        with a dict of atmos and recon (or just atmos).  The question is how to handle
        this for the various children...
        Thats why the idstr can be passed.  The idstr is used to find the correct resource sharing object, and
        then it is given the parents.  Note, that this module will not allow 2 resource sharing objects to have
        the same idstr.
        """
        if idstr not in self.idstr:
            raise Exception("xinterp_dm - newParent() idstr not known. %s %s"%(idstr,self.idstr))
        indx=self.idstr.index(idstr)
        if type(parent)!=type({}):
            parent={"1":parent}
        self.thisObjList[indx].parent=parent
        
    def generateNext(self,ms=None):
        """DM main loop:  Reflect from DM (add DM figure to phase)"""
        this=self.thisObjList[self.currentIdObjCnt]
        if not this.parent.has_key("atmos") and not this.parent.has_key("recon"):
            print "zdm object assigning parents automatically"
            keylist=this.parent.keys()
            for key in keylist:
                if this.parent[key].outputData.shape[-2:]==(self.npup,self.npup):#the right shape for an atmos object
                    print "zdm parent object %s becoming atmos with output shape %s"%(str(key),str(this.parent[key].outputData.shape))
                    this.parent["atmos"]=this.parent[key]
                else:
                    print "zdm parent object %s becoming recon with output shape %s"%(str(key),str(this.parent[key].outputData.shape))
                    this.parent["recon"]=this.parent[key]
            for key in keylist:
                del(this.parent[key])
            if not this.parent.has_key("atmos"):#maybe just running for poking?
                print "zdm object no parent atmos object found.  Assuming unperturbed phase in"

        if self.generate==1:
            if self.newDataWaiting:
                if this.parent.has_key("atmos"):
                    if this.parent["atmos"].dataValid==1:
                        self.outputData[:]=this.parent["atmos"].outputData
                        self.dataValid=1
                        if self.control["phaseCovariance"]:
                            self.montePhaseCovariance()
                        elif self.lastPhaseCovariance:
                            self.finishMontePhaseCovariance()
                        self.lastPhaseCovariance=self.control["phaseCovariance"]
                    else:
                        print "DM: waiting for data form atmos, but not valid"
                        self.dataValid=0
                else:#no atmos parent
                    self.dataValid=1#zero output always valid.
                if self.control["dm_update"]==1 and self.currentIdObjCnt==0 and this.parent["recon"].dataValid==1:
                    if self.actStart==None:
                        self.getActuatorOffsets(this.parent["recon"].outputData)
                    self.reconData=this.parent["recon"].outputData[self.actStart:self.actEnd]
                    if self.dmInfo.iirCoeffs is not None:
                        self.applyIIR(self.reconData)
                    self.update()
                    self.dataValid=1#update the output.
            if self.dataValid:
                if self.sendFullDM:#output full surface
                    self.outputData=self.dmphs
                else:
                    if this.parent.has_key("atmos"):
                        addToOutput=1
                    else:
                        addToOutput=0
                    this.lineOfSight.selectSubPupil(self.outputData,addToOutput,removePiston=1)
                    self.selectedDmPhs=this.lineOfSight.selectedDmPhs

                """
                if self.subpxlInterp:
                    #do the interpolation...
                    if this.xoffsub==0 and this.yoffsub==0:#no interp needed
                        pass
                    else:
                        gslCubSplineInterp(this.selectedDmPhs,self.yaxisInterp,self.xaxisInterp,
                                           this.yaxisInterp,this.xaxisInterp,self.interpolated,
                                           self.interpolationNthreads)
                        self.selectedDmPhs=self.interpolated

                if this.wavelengthAdjustor==1:#dm is shaped for this wavelength...
                    if this.parent.has_key("atmos"):
                        self.outputData+=self.selectedDmPhs
                    else:
                        self.outputData[:,]=self.selectedDmPhs
                else:
                    if this.parent.has_key("atmos"):
                        self.outputData+=self.selectedDmPhs*this.wavelengthAdjustor
                    else:
                        self.outputData[:,]=self.selectedDmPhs*this.wavelengthAdjustor
                        
                if self.pupil!=None:
                    self.outputData*=self.pupil.fn
                    #now remove piston
                    phasesum=numpy.sum(self.outputData.ravel())
                    self.outputData-=phasesum/self.pupil.sum
                    self.outputData*=self.pupil.fn
                #self.outputData += (self.dmphs*(self.wavelengthRatio))# Calculate reflected phase
                """
        else:
            self.dataValid=0

    def applyIIR(self,reconData):
        """Applies an IIR filter to the reconstructor data"""
        a=self.dmInfo.iirCoeffs[0]
        b=self.dmInfo.iirCoeffs[1]
        for i in range(len(reconData)):
            new_mode_val, self.dmInfo.iirZi[i] = lfilter(b, a, [reconData[i]], zi=self.dmInfo.iirZi[i])
            self.reconData[i] = new_mode_val[0]
        #
        # print self.reconData
        
    def getActuatorOffsets(self,reconOutput):
        if self.nmodes==reconOutput.size:
            self.actStart=0
            self.actEnd=reconOutput.size
        else:#find out which part of the recon output we should be using.
            self.dmList=self.dmObj.makeDMList(self.dmInfo.actuatorsFrom)
            nactsList=[]
            nactsCumList=[0]
            for dm in self.dmList:
                if dm.zonalDM:
                    if self.pupil==None:
                        r2=0
                    else:
                        r2=self.pupil.r2
                    tmp=dm.computeDMPupil(self.atmosGeom,centObscuration=r2,retPupil=0)
                    # tmp is dmflag,subarea (or None,None for modal DMs.)
                    #self.dmPupList.append(tmp[0])

                    nactsList.append(int(numpy.sum(tmp[0].ravel())))
                    #if dm.pokeSpacing!=None:
                    #    self.npokesList.append(dm.pokeSpacing**2)
                    #else:
                    #    self.npokesList.append(self.nactsList[-1])
                    #self.npokesCumList.append(self.npokesCumList[-1]+self.npokesList[-1])
                else:#a modal DM
                    #self.dmPupList.append(None)
                    nactsList.append(dm.nact)#nact is the number of modes
                    #self.npokesList.append(dm.nact)
                    #self.npokesCumList.append(self.npokesCumList[-1]+dm.nact)
                nactsCumList.append(nactsList[-1]+nactsCumList[-1])
                #self.closedLoopList.append(dm.closedLoop)
            print self.dmList
            print self.dmInfo
            indx=self.dmList.index(self.dmInfo)
            self.actStart=nactsCumList[indx]
            self.actEnd=nactsCumList[indx+1]
            print "nactsCumList computed by DM: %s"%nactsCumList
        print "DM %s using actuators from %d to %d (recon output size %d)"%(self.dmInfo.label,self.actStart,self.actEnd,reconOutput.size)

            
    def update(self):
        """DM figure update: Add zernike updates to
        mirror (if the loop is closed)"""
	self.dmphs[:]=0#numpy.zeros((self.npup,self.npup),numpy.float64)
        #zoffset is only used here (I think).  Should be set by GUI.
        #zpoke also used by reconstructor (but doesn't do anything there)
        #powers=self.gamma*self.reconData
	for j in range(1,self.nmodes):
	    self.dmphs+=self.reconData[j]*self.zern[j,:,:]     # add closed loop corrections
	    
        if type(self.zoffset)!=type(None):
            for j in range(min(self.nmodes,self.zoffset.shape[0])):
                self.dmphs+=self.zoffset[j]*self.zern[j,:,:]#always add constant offsets or pokes
                

    def montePhaseCovariance(self):
        a=numpy.zeros((self.nmodes,),numpy.float64)
        phs=self.parent["atmos"].outputData
        #if type(self.mirrorModes)==type(None):
        #    self.makeMirrorModes(fitpup=self.fitPupMirrorModes)
        start=1
        for i in xrange(start,self.nmodes):
            a[i]=numpy.sum(phs*self.zern[i])
        if type(self.monteNoll)==type(None):
            self.monteNoll=numpy.zeros((self.nmodes,self.nmodes),numpy.float64)
            self.navNoll=0
        self.monteNoll+=a[:,numpy.newaxis]*a[numpy.newaxis,:]
        self.navNoll+=1
    def finishMontePhaseCovariance(self):
        self.montePhaseCovMat=self.monteNoll/(self.navNoll*self.zern[0].sum()**2)#*self.pupil.sum)#*self.pxlscale**2)#should pupil.sum be with or without the central obscuration? (not that it makes a huge difference!).  Don't scale to pixels, because the poke matrix scales from radians to pixels.
        if self.montePhaseCovfname!=None:
            print "%s has been written"%self.montePhaseCovfname
            util.FITS.Write(self.montePhaseCovMat,self.montePhaseCovfname)
        #This should be similar to self.Noll (check scaling - am not sure its quite right).


    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        this=self.thisObjList[self.sentPlotsCnt]
        if this.idstr==None or this.idstr=="":
            id=""
        else:
            id=" (%s)"%this.idstr
        txt=""
        if self.sentPlotsCnt==0:
            txt+="""<plot title="zdm output%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="zdm mirror%s" cmd="data=-%s.dmphs" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt += """<plot title="zdm modes%s" cmd="data=%s.reconData" ret="data" type="pylab" when="rpt"/>\n""" % (
            id, objname)
        txt+="""<plot title="zdm selected mirror%s" cmd="data=-%s.thisObjList[%d].lineOfSight.selectedDmPhs" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt)
        self.sentPlotsCnt=(self.sentPlotsCnt+1)%len(self.thisObjList)
        return txt

    def getParams(self):
        """parameters required for this module, in the form of {"paramName":defaultValue,...}
        These params can then be placed in the config file... if not set by the
        user, the param should still be in config file as default value for
        future reference purposes.
        """
        #This is a working example.  Please feel free to change the parameters
        #required. (if you do, also change the config.getParam() calls too).
        paramList=[]
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))
        paramList.append(base.dataType.dataType(description="dmObj",typ="code",val="import util.dm;dmObj=util.dm.dmOverview(dmInfoList,atmosGeom=atmosGeom)",comment="TODO: dmInfoList is a list of util.dm.dmInfo objects, initialised with (label (for this dm), idlist (list of (dm ID,source ID) or just sourceID, the idstr for a particular DM object (height and direction) and the idstr for a given source direction), height, nact, and various other things for which the defaults may be okay (see the source code))"))
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(npup,npup/2,npup/2*telSec/telDiam,wfs_nsubx,wfs_minarea)",comment="Telescope pupil object"))
        paramList.append(base.dataType.dataType(description="zoffset",typ="eval",val="None",comment="Arbitrary zernike offset"))






### Startup ############################################################
if __name__=="__main__":
    print "Not yet implemented correctly..."
    raise Exception("Debug running not yet implemented correctly")




