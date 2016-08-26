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
#from cmod.interp import gslCubSplineInterp
import base.aobase
import util.FITS,util.dist
import util.zernikeMod
from util.dm import MirrorSurface

class dm(base.aobase.aobase):
    """
    Xinetics DM simulation object:
    Bicubically or spline interpolated actuator map.
    Add DM phases to input phase array and pass to output
    Update mirror figure on request from reconstructor (when data becomes
    available)
    This object is able to implement resource sharing, but only to the extent
    of having one reconstructor parent - ie it is able to physically represent
    one DM, taking atmos from lots of sky directions.  
    The mirror shape itself is only updated for the first resource sharing
    object, and then stays like this for all the others.
    Note that this DM object can only have 1 unique conjugate height -
    otherwise it would be unphysical...

    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        if type(parent)!=type({}):
            parent={"1":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.datatype=self.config.getVal("xinterpdmDataType",default=numpy.float32)
        #self.sendFullDM=self.config.getVal("sendFullDM",default=0)#used if connecting to wideField.py science module
        if forGUISetup==1:
            #if self.sendFullDM:
            self.dmObj=self.config.getVal("dmOverview",default=self.config.getVal("dmObj"),raiseerror=0)
            if self.dmObj.getDM(idstr).sendFullDM:
                dmpup=self.dmObj.calcdmpup(self.idstr[0])
                self.outputData=[(dmpup,dmpup),"f"]
            else:
                npup=self.config.getVal("npup")
                self.outputData=[(npup,npup),self.datatype]
        else: # set up for simulation.
            self.control={"dm_update":1,"zoffset":None,"phaseCovariance":0}#,"zpoke":numpy.zeros((self.nact*self.nact,),self.datatype)}#,"poke":0}
            
            # Extra data for interpolated DM 
            #self.actmap_1d=numpy.zeros((self.nact*self.nact),self.datatype)
            #self.pokemap=numpy.zeros((self.nact*self.nact),self.datatype)
            #self.pokeval=self.config.getVal("pokeval")

            self.atmosGeom=self.config.getVal("atmosGeom",default=None,raiseerror=0)
            self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
            if self.dmObj==None or type(self.dmObj)!=type(self.atmosGeom):
                print "DEPRECATION: warning: dmObj should now be dmOverview"
                self.dmObj=self.config.getVal("dmObj",default=None,raiseerror=0)
            self.npup=self.config.getVal("npup")
            self.pupil=self.config.getVal("pupil",raiseerror=0)
            self.interpolationNthreads=self.config.getVal("interpolationNthreads", default=0)
            self.subtractTipTilt=self.config.getVal("subtractTipTilt",default=0)
            self.monteNoll=None
            self.mirrorModes=None
            self.mirrorScale=None
            self.actStart=None#offsets into the parent actutator lists...
            self.actEnd=None
            self.navNoll=0
            self.allZero=0
            self.lastPhaseCovariance=0
            if self.dmObj==None:#depreciated mode
                self.rotation=None
                self.dmTiltAngle=0.
                self.dmTiltTheta=0.
                self.subpxlInterp=self.config.getVal("subpxlInterp",default=0)#should we do sub pixel interpolation for non-zero conjugated DMs?
                self.nact=self.config.getVal("nAct")
                self.dmpup=self.config.getVal("dmpup",default=self.npup)#number of pixels to use to store the DM phase - this may be larger than npup if the DM is conjugated somewhere other than ground. depreciated
                self.conjHeight=self.config.getVal("dmConjugateHeight",default=0.)#the height at which this pupil phase should be conjugated too. depreciated

                print "xinterp_dm - dmObj not found in param file - continuing in depreciated mode"
                sourceID=self.config.getVal("sourceID",default=self.idstr[0])
                if self.atmosGeom==None:
                    self.sourceLam=self.config.getVal("sourceLam")
                else:
                    self.sourceLam=self.atmosGeom.sourceLambda(sourceID)
                    if self.sourceLam==None:
                        self.sourceLam=self.config.getVal("sourceLam")
                        print "WARNING - sourceLam not found in atmosGeom, using %g"%self.sourceLam
                self.dmflag=self.config.getVal("dmflag")#depreciated
                self.actCoupling=config.getVal("actCoupling")#coupling of actuators to neighbours, eg, 0.1. depreciated
                self.interpType=self.config.getVal("dmInterpType",default="spline")
                self.actSlaves=None
                self.stuckActs=config.getVal("stuckActs",default=None,raiseerror=0)
                self.actFlattening=config.getVal("actFlattening",default=1.)#flattening of gradients, default 1., should be between 0-1.
                self.mirrorSurface=MirrorSurface(self.interpType,self.dmpup,self.nact,1,
                                                 self.actoffset,self.actCoupling,self.actFlattening,
                                                 interpolationNthreads = self.interpolationNthreads,stuckActs=self.stuckActs)
                self.maxStroke=0#depreciated mode dones't have max stroke.
                self.dmDynamics=None
            else:
                self.thisdm=self.dmObj.getDM(self.idstr[0])
                self.dmTiltAngle=self.thisdm.tiltAngle
                self.dmTiltTheta=self.thisdm.tiltTheta
                self.rotation=self.thisdm.rotation#rotation can be an angle to rotate the DM by, or a function which when evaluated (every iteration) returns a new angle - eg to use for a rotating mirror.

                self.interpType=self.thisdm.interpType
                self.alignmentOffset=self.thisdm.alignmentOffset#to simulate error in dm alignment, tuple of x,y pixels offset.
                self.subpxlInterp=self.thisdm.subpxlInterp
                self.dmpup=self.dmObj.calcdmpup(self.idstr[0])
                self.conjHeight=self.dmObj.getHeight(self.idstr[0])
                self.actoffset=self.dmObj.getactoffset(self.idstr[0])
                self.nact=self.dmObj.getnAct(self.idstr[0])
                self.sourceLam=self.thisdm.reconLam#the reconstructor wavelength - ie wavelength for which the DM is shaped.
                #thisdm.maxStroke is in microns.  So, convert to radians, and divide by 2, so that have half in each direction.
                self.maxStroke=self.thisdm.maxStroke*1000./self.sourceLam*2*numpy.pi/2.
                #sourceID=self.dmObj.getSourceID(self.idstr[0])
                # need a flag telling us which actuators are valid to use (otherwise we waste time and memory reconstructing all actuators).  This flag will depend on the geometry (Freid etc) and dmpup.
                if self.pupil is None:
                    r2=0
                else:
                    r2=self.pupil.r2
                self.dmflag=self.dmObj.computeDMPupil(self.idstr[0],centObscuration=r2,retPupil=0)[0]
                self.actCoupling=self.dmObj.getcoupling(self.idstr[0])
                self.actSlaves=self.thisdm.getSlaving()
                self.mirrorSurface = self.thisdm.getMirrorSurface(phsOut = 1,                                                                interpolationNthreads = self.interpolationNthreads)
                self.dmDynamics=self.thisdm.dmDynamics#an array of the fraction of shift to new position that occur each timestep, e.g. for a simulation with the WFS updating every 4 frames, this could be [0.5,0.5,0.5,1.] would move 50% after 1 step, 75% after 2 steps, 87.5% after 3 steps, and arrive after 4 steps.
            self.lastactmap=None#only used if dmDynamics are in use.
            self.dynamicStep=0#only used if dmDynamics are in use.
            if self.subtractTipTilt:
                self.tilt=numpy.arange(self.nact)-self.nact/2.+0.5
                self.tilt/=numpy.sqrt((self.tilt**2).sum()*self.nact)
            if self.conjHeight==0 and self.dmTiltAngle==0:
                self.subpxlInterp=0#no need for interpolation...
                
            if self.subpxlInterp or self.alignmentOffset[0]!=0 or self.alignmentOffset[1]!=0:
                #self.interpolated=numpy.zeros((self.npup,self.npup),numpy.float32)
                self.yaxisInterp=numpy.arange(self.npup+1).astype(numpy.float64)#Must be float 64 because of gsl restriction (spline functions require it)
                self.xaxisInterp=numpy.arange(self.npup+1).astype(numpy.float64)
            else:
                self.yaxisInterp=None
                self.xaxisInterp=None
                #self.interpolated=None
            self.actmap=numpy.zeros((self.nact,self.nact),self.datatype)
            #self.nsubx=n =self.config.getVal("wfs_nsubx")
            #self.wfsn=self.config.getVal("wfs_n")
            #self.dmminarea=self.config.getVal("dmminarea",default=0.25)
            self.dmphs=numpy.zeros((self.dmpup,self.dmpup),numpy.float32) # DM figure
            self.mirrorSurface.phsOut=self.dmphs
            #we now need to work out which part of the DM to use... this depends on conjugate height, and source location.
            self.reconData=None#numpy.zeros(self.nact*self.nact,self.datatype)
            

            self.nacts=int(numpy.sum(numpy.sum(self.dmflag)))
            print "xinterp_dm_%s: using nacts=%d"%(str(self.idstr[0]),self.nacts)
            self.fitPupMirrorModes=self.config.getVal("fitPupMirrorModes",default=1)
            #self.subflag=self.config.getVal("subflag")
            #self.dmflag_1d=numpy.reshape(self.dmflag,(self.nact*self.nact,))
            #self.subflag_1d=numpy.reshape(self.subflag,(n*n,))
            #self.wfsdata = numpy.sum(self.subflag_1d) #Number of filled subaps
            #self.dmdata=numpy.sum(self.dmflag_1d) # Number of used actuators
            self.dmindices=numpy.nonzero(self.dmflag.ravel())[0]
            self.telDiam=self.config.getVal("telDiam")
            self.sendFullDM=self.dmObj.getDM(idstr).sendFullDM
            
            if self.sendFullDM:
                self.outputData=self.dmphs
            else:
                self.outputData=numpy.zeros((self.npup,self.npup),self.datatype)
            self.lowOrderModeDict={}#dict containing low order modes which can be subtracted from the atmos phase.
            self.lowOrderModeNormDict={}
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
        
        this.lineOfSight=util.dm.DMLineOfSight(self.dmpup,self.npup,self.conjHeight,self.dmphs,this.sourceAlt,this.sourceTheta,this.sourcePhi,self.telDiam,wavelengthAdjustor,self.xaxisInterp,self.yaxisInterp,self.dmTiltAngle,self.dmTiltTheta,self.alignmentOffset,self.subpxlInterp,self.pupil,self.interpolationNthreads)

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
        if self.atmosGeom==None:
            this.sourceAlt=this.config.getVal("sourceAlt")#depreciated
            this.sourceLam=this.config.getVal("sourceLam")          # Wavelength for this optical path
        else:
            if self.dmObj==None:
                this.sourceID=this.config.getVal("sourceID",default=idstr)#sourceID is the idstr for the source direction as provided for the atmosGeom object. depreciated
            else:
                this.sourceID=self.dmObj.getSourceID(idstr)
            this.sourceAlt=self.atmosGeom.sourceAlt(this.sourceID)
            this.sourceLam=self.atmosGeom.sourceLambda(this.sourceID)#this.config.getVal("sourceLam")          # Wavelength for this optical path
            if this.sourceLam==None:
                this.sourceLam=this.config.getVal("sourceLam")          # Wavelength for this optical path
                print "WARNING - sourceLam not found in atmosGeom, using %g"%this.sourceLam

        this.xaxisInterp=None
        if self.conjHeight!=0 or self.dmTiltAngle!=0:
            if this.sourceAlt>0 and self.subpxlInterp==0:
                raise Exception("xinterp_dm - finite altitude source specified but subpxlInterp==0.")
            if self.atmosGeom==None:
                this.sourceTheta=this.config.getVal("sourceTheta")*numpy.pi/180/3600.#radians depreciated
                this.sourcePhi=this.config.getVal("sourcePhi")*numpy.pi/180#radians depreciated
            else:
                this.sourceTheta=self.atmosGeom.sourceTheta(this.sourceID)*numpy.pi/180/3600.
                this.sourcePhi=self.atmosGeom.sourcePhi(this.sourceID)*numpy.pi/180
            r=self.conjHeight*numpy.tan(this.sourceTheta)/self.telDiam*self.npup
            this.xoff+=self.alignmentOffset[0]
            this.xoffend+=self.alignmentOffset[0]
            this.yoff+=self.alignmentOffset[1]
            this.yoffend+=self.alignmentOffset[1]
            if self.subpxlInterp:
                if this.sourceAlt>0:
                    if this.sourceAlt<self.conjHeight:
                        #source is below mirror height...???
                        
                        #I think this just means that we need to flip the mirror actuators...
                        print "Warning - you have a DM conjugated above a source height - is this what you expect? (actuators will be reversed!)"
                    #Compute factors for a tilted DM...
                    xfact=numpy.cos(self.dmTiltTheta*numpy.pi/180.)*(1/numpy.cos(self.dmTiltAngle*numpy.pi/180.)-1)+1
                    yfact=numpy.sin(self.dmTiltTheta*numpy.pi/180.)*(1/numpy.cos(self.dmTiltAngle*numpy.pi/180.)-1)+1

                    x=r*numpy.cos(this.sourcePhi)
                    y=r*numpy.sin(this.sourcePhi)
                    #xoff+x would be the starting point if source at infinity.
                    #xoff+x+npup would be the ending point.
                    #If a tilted mirror, xoff+x-(npup*xfact-npup)/2 would be starting point.
                    #So, I need to reduce this by a factor of (alt-conj)/alt
                    npxl=self.npup*(this.sourceAlt-self.conjHeight)/this.sourceAlt
                    #if npxl <0, means have to flip the mirror...
                    #this.xoff=(this.xoff+x)+(npup-numpy.fabs(npxl))/2.-(npxl*xfact-npup)/2.
                    #this.yoff=(this.yoff+y)+(npup-numpy.fabs(npxl))/2.-(npxl*yfact-npup)/2.
                    this.xoff=(this.xoff+x)+(npup-numpy.fabs(npxl*xfact))/2.
                    this.yoff=(this.yoff+y)+(npup-numpy.fabs(npxl*yfact))/2.
                    this.xoffsub=this.xoff%1
                    this.yoffsub=this.yoff%1
                    tol=1e-10#allow for floating point inprecision 
                    this.xoffend=int(numpy.ceil(this.xoff+numpy.fabs(npxl)*xfact-tol))#oversize by 1pxl
                    this.yoffend=int(numpy.ceil(this.yoff+numpy.fabs(npxl)*yfact-tol))#unless exact fit
                    this.xoff=int(this.xoff+tol)
                    this.yoff=int(this.yoff+tol)
                    if npxl>0:#source above DM conjugate...
                        this.xaxisInterp=numpy.arange(self.npup).astype(numpy.float64)*(npxl-1.)/(self.npup-1.)*xfact+this.xoffsub
                        this.yaxisInterp=numpy.arange(self.npup).astype(numpy.float64)*(npxl-1.)/(self.npup-1.)*yfact+this.yoffsub
                    else:#need to flip since will see mirror in reverse (I think)
                        this.xaxisInterp=numpy.arange(self.npup-1,-1,-1).astype(numpy.float64)*(numpy.fabs(npxl)-1.)/(self.npup-1.)*xfact+this.xoffsub
                        this.yaxisInterp=numpy.arange(self.npup-1,-1,-1).astype(numpy.float64)*(numpy.fabs(npxl)-1.)/(self.npup-1.)*yfact+this.yoffsub
                else:#source at infinity.
                    this.xoff+=r*numpy.cos(this.sourcePhi)
                    this.yoff+=r*numpy.sin(this.sourcePhi)
                    if self.dmTiltAngle==0:
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
                    else:#DM is tilted... this will change effective actuator spacing, and hence the part of the DM that we interpolate...
                        #Instead of seeing N actuators, you will see nN actuators where n>1, ie an increased number of actuators.
                        #The increased effective size is oldsize/cos(theta), ie this is the number of actuators that will be seen
                        #So the additional number is oldsize(1/cos(theta)-1)
                        #So, nact*((1/cos(theta)-1)*cosorsin(theta)+1) is the new number of actuators...
                        xfact=numpy.cos(self.dmTiltTheta*numpy.pi/180.)*(1/numpy.cos(self.dmTiltAngle*numpy.pi/180.)-1)+1
                        yfact=numpy.sin(self.dmTiltTheta*numpy.pi/180.)*(1/numpy.cos(self.dmTiltAngle*numpy.pi/180.)-1)+1
                        #shift the start backwards.
                        this.xoff-=(npup*xfact-npup)/2.
                        this.xoffsub=this.xoff%1
                        this.xoffend=int(numpy.ceil(this.xoff+self.npup*xfact))
                        this.xoff=int(this.xoff)
                        this.xaxisInterp=numpy.arange(self.npup).astype(numpy.float64)*xfact+this.xoffsub
                        this.yoff-=(npup*yfact-npup)/2.
                        this.yoffsub=this.yoff%1
                        this.yoff=int(this.yoff)
                        this.yoffend=this.yoff+int(numpy.ceil(this.yoffsub+self.npup*yfact))
                        this.yaxisInterp=numpy.arange(self.npup).astype(numpy.float64)*yfact+this.yoffsub
                #Need to check that self.xaxisInterp is big enough...
                if self.yaxisInterp.shape[0]<(this.yoffend-this.yoff):
                    self.yaxisInterp=numpy.arange(this.yoffend-this.yoff).astype(numpy.float64)
                    print "Increasing yaxisInterp"
                if self.xaxisInterp.shape[0]<(this.xoffend-this.xoff):
                    self.xaxisInterp=numpy.arange(this.xoffend-this.xoff).astype(numpy.float64)
                    print "Increasing xaxisInterp"
            else:#not doing sub pixel interpolation
                if self.dmTiltTheta!=0:
                    raise Exception("not doing subpixel interpolation but tilt specified for DM")
                this.xoff+=int(r*numpy.cos(this.sourcePhi)+0.5)#round
                this.yoff+=int(r*numpy.sin(this.sourcePhi)+0.5)#correctly
                this.xoffend=this.xoff+self.npup
                this.yoffend=this.yoff+self.npup
        if this.xaxisInterp==None:#ground conjugated...
            if self.alignmentOffset[0]!=0 or self.alignmentOffset[1]!=0:
                this.xaxisInterp=numpy.arange(self.npup).astype(numpy.float64)+self.alignmentOffset[0]
                this.yaxisInterp=numpy.arange(self.npup).astype(numpy.float64)+self.alignmentOffset[1]
            #but don't bother with any interpolation here...
        # select the correct part of the dm phs.
        if this.xoff<0 or this.yoff<0 or this.xoffend>self.dmpup or this.yoffend>self.dmpup:
            raise Exception("DM pupil not large enough to hold all the conjugated sources... %d %d %d %d %d"%(this.xoff,this.yoff,this.xoffend,this.yoffend,self.dmpup))

        this.selectedDmPhs=self.dmphs[this.yoff:this.yoffend,this.xoff:this.xoffend]
        #this.wfsLam=this.config.getVal("wfslam",default=this.sourceLam)
        #this.wavelengthRatio=this.wfsLam/this.sourceLam#wfs_lam/lam
        this.wavelengthAdjustor=self.sourceLam/this.sourceLam#this.wavelengthRatio/self.wavelengthRatio#the mirror will be shaped as for self.wavelengthRatio...if this.sourceLam is longer than self.sourceLam, radians of phase P-V will be smaller, so less change needed in wavelengthAdjustor.

        #self.tempphs=numpy.zeros((self.npup,self.npup),self.datatype,savespace=1)            # DM figure      
        #self.phs=numpy.zeros((self.npup,self.npup),self.datatype)                # Output phase array
        #self.tilt_gain=config.getVal("tilt_gain")
##         self.gamma=config.getVal("gamma")
##         self.useVariableGamma=config.getVal("variableGamma",default=0)
##         if self.useVariableGamma:
##             gamma=numpy.zeros((self.nact,self.nact),numpy.float32)
##             for i in range(self.nact):
##                 for j in range(self.nact):
##                     xs=j*self.wfsn-self.wfsn/2
##                     xe=xs+self.wfsn
##                     ys=i*self.wfsn-self.wfsn/2
##                     ye=ys+self.wfsn
##                     if xs<0:xs=0
##                     if ys<0:ys=0
##                     if xe>self.npup:xe=self.npup
##                     if ye>self.npup:ye=self.npup
##                     s=numpy.sum(numpy.sum(self.pupil.fn[ys:ye,xs:xe]))/float(self.wfsn*self.wfsn)
##                     gamma[i,j]=self.gamma*s**0.
##             self.gamma=gamma
        """
    def finalInitialisation(self):
        """Just check that there aren't 2 objects with same idstr..."""
        tmp=[]
        for id in self.idstr:
            if id in tmp:
                raise Exception("xinterp_dm - resource sharing cannot have 2 objects with same idstr")
            tmp.append(id)
        
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

            
### DM main loop ########################################################

# reflect from DM (add DM figure to phs)

        
    def generateNext(self,ms=None):
        """DM main loop:  Reflect from DM (add DM figure to phase)"""

        this=self.thisObjList[self.currentIdObjCnt]
        if not this.parent.has_key("atmos") and not this.parent.has_key("recon"):
            print "xinterp_dm object assigning parents automatically"
            keylist=this.parent.keys()
            for key in keylist:
                if this.parent[key].outputData.shape==(self.npup,self.npup):#the right shape for an atmos object
                    print "xinterp_dm parent object %s becoming atmos"%str(key)
                    this.parent["atmos"]=this.parent[key]
                else:
                    print "xinterp_dm parent object %s becoming recon with output shape %s"%(\
                        str(key),str(this.parent[key].outputData.shape))
                    this.parent["recon"]=this.parent[key]
            for key in keylist:
                del(this.parent[key])
            if not this.parent.has_key("atmos"):#maybe we're just using this for a poke matrix?
                print "xinterp_dm object no parent atmos object found.  Assuming unperturbed phase in."
        if self.generate==1:
            if self.newDataWaiting:
                if this.parent.has_key("atmos"):
                    if this.parent["atmos"].dataValid==1:
                        self.outputData[:,]=this.parent["atmos"].outputData # make a copy
                        self.dataValid=1
                        if self.control["phaseCovariance"]:
                            self.montePhaseCovariance()
                        elif self.lastPhaseCovariance:
                            self.finishMontePhaseCovariance()
                        self.lastPhaseCovariance=self.control["phaseCovariance"]
                        # now remove any low order modes as requested from the atmos - eg tip.tilt for lgs
                        #for mode in this.subLowOrderModeList:
                        #    coeff=numpy.sum(numpy.sum(self.outputData*self.lowOrderModeDict[mode]/self.lowOrderModeNormDict[mode]))
                        #    self.outputData-=coeff*self.lowOrderModeDict[mode]
                    else:
                        print "DM: waiting for data from atmos, but not valid"
                        self.dataValid=0
                else:#has not atmos parent...
                    self.dataValid=1#zero output always valid...
                if self.control["dm_update"]==1 and self.currentIdObjCnt==0:
                    if this.parent["recon"].dataValid==1:
                        #we are the first resource sharer, and have reconstructor commands ready...
                        if self.actStart==None:
                            self.getActuatorOffsets(this.parent["recon"].outputData)
                        self.reconData=this.parent["recon"].outputData[self.actStart:self.actEnd]
                        self.dynamicStep=0
                        self.update()    # update the dm figure.
                        self.dataValid=1 # update the output.
                    elif self.dmDynamics!=None:
                        self.dynamicStep+=1
                        self.update()
                        self.dataValid=1
            if self.dataValid:
                if self.sendFullDM:#output full DM surface, not just along one line of sight.  Used by wideField.py science module
                    self.outputData=self.dmphs
                elif not self.allZero:
                    if this.parent.has_key("atmos"):
                        addToOutput=1
                    else:
                        addToOutput=0
                    this.lineOfSight.selectSubPupil(self.outputData,addToOutput,removePiston=1)
                    self.selectedDmPhs=this.lineOfSight.selectedDmPhs

                    """
                    self.selectedDmPhs=this.selectedDmPhs
                    if self.subpxlInterp:
                        #do the interpolation...
                        if this.xoffsub==0 and this.yoffsub==0:#no interp needed
                            pass
                        else:
                            gslCubSplineInterp(this.selectedDmPhs,self.yaxisInterp,self.xaxisInterp,
                                                  this.yaxisInterp,this.xaxisInterp,self.interpolated,
                                                  self.interpolationNthreads)
                            self.selectedDmPhs=self.interpolated
                    elif self.alignmentOffset[0]!=0 or self.alignmentOffset[1]!=0:
                        #ground conjugate or not interpolating, but we want to offset anyway...
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
                        #now remove any piston... though probably not necessary.
                        phasesum=numpy.sum(self.outputData.ravel())
                        self.outputData-=phasesum/self.pupil.sum
                        self.outputData*=self.pupil.fn
                    """
                else:
                    if not this.parent.has_key("atmos"):
                        self.outputData[:]=0
                    #else:
                    #    this.lineOfSight.selectSubPupil(self.outputData,1,removePiston=1)
                    #    self.selectedDmPhs=this.lineOfSight.selectedDmPhs
        else:
            self.dataValid=0
            
    def getActuatorOffsets(self,reconOutput):
        if self.nacts==reconOutput.size:
            self.actStart=0
            self.actEnd=reconOutput.size
        else:#find out which part of the recon output we should be using.
            self.dmList=self.dmObj.makeDMList(self.thisdm.actuatorsFrom)
            nactsList=[]
            nactsCumList=[0]
            for dm in self.dmList:
                if dm.zonalDM:
                    if self.pupil==None:
                        r2=0
                    else:
                        r2=self.pupil.r2
                    tmp=dm.getDMFlag(self.atmosGeom,centObscuration=r2)
                    nactsList.append(int(tmp.sum()))
                else:#a modal DM
                    #self.dmPupList.append(None)
                    nactsList.append(dm.nact)#nact is the number of modes
                    #self.npokesList.append(dm.nact)
                    #self.npokesCumList.append(self.npokesCumList[-1]+dm.nact)
                nactsCumList.append(nactsList[-1]+nactsCumList[-1])
                #self.closedLoopList.append(dm.closedLoop)
            print self.dmList
            print self.thisdm
            indx=self.dmList.index(self.thisdm)
            self.actStart=nactsCumList[indx]
            self.actEnd=nactsCumList[indx+1]
            print "nactsCumList computed by DM: %s"%nactsCumList
        print "DM %s using actuators from %d to %d (recon output size %d)"%(self.thisdm.label,self.actStart,self.actEnd,reconOutput.size)

    def update(self):
        """DM figure update: Add actuator updates to
        mirror (if the loop is closed)"""

        #zoffset is only used here (I think).  Should be set by GUI.
        zoffset=self.control["zoffset"]#self.controlDict['zoffset']
        self.allZero=0
        if type(zoffset)==type(None):
            if numpy.alltrue(self.reconData==0):
                if self.allZero==0:
                    self.dmphs[:]=0
                self.allZero=1#this saves some computation time when poking tomographic systems...
            self.actmap[:,]=0
        else:
            self.actmap[:,]=numpy.reshape(zoffset,(self.nact,self.nact))
        #zpoke also used by reconstructor (but doesn't do anything there)
        #zpoke=self.control["zpoke"]#self.controlDict['zpoke']
        if not self.allZero:
            l=len(self.reconData.shape)
            if l==3:#output probably from SOR module...
                self.actmap+=self.reconData[0]
                self.geom="hudgin"
            elif l==2:
                self.actmap+=self.reconData
                self.geom="fried"
            elif l==1:#output from the xinterp_recon/tomoRecon module...
                if self.reconData.shape[0]==self.nact*self.nact:
                    self.actmap+=numpy.reshape(self.reconData,(self.nact,self.nact))
                else:
                    numpy.put(self.actmap.ravel(),self.dmindices,self.reconData)
                    #can do actuator slaving here.
                    if self.actSlaves!=None:
                        self.applySlaving(self.actmap.ravel(),self.actSlaves)
                self.geom=None#"fried"#use the actoffset instead.  If this is zero (default), same as fried.
            if self.dmDynamics!=None:
                #actmap has the shape we want to move to.  lastactmap has the current shape.  self.dynamicStep has the step number.
                frac=self.dmDynamics[self.dynamicStep]
                if self.lastactmap==None:
                    self.lastactmap=numpy.zeros(self.actmap.shape,self.actmap.dtype)
                self.actmap[:]=self.lastactmap+frac*(self.actmap-self.lastactmap)
            if self.subtractTipTilt:
                #remove tip and tilt from actmap.
                #Do this in a lazy way... which isn't totally accurate for non-square geometry (i.e. ie uses corner actuators too).
                t=(self.actmap.sum(0)*self.tilt).sum()
                self.actmap-=t*self.tilt
                t=(self.actmap.sum(1)*self.tilt).sum()
                self.actmap.T-=t*self.tilt
            if self.maxStroke!=0:#maxStrke is in radians at sourceLam wavelength.
                #remove mean
                self.actmap-=self.actmap.mean()
                self.actmap[:]=numpy.where(self.actmap>self.maxStroke,self.maxStroke,self.actmap)
                self.actmap[:]=numpy.where(self.actmap<-self.maxStroke,-self.maxStroke,self.actmap)
            self.mirrorSurface.fit(self.actmap)
            if self.dmDynamics!=None:#save the state for next time
                self.lastactmap[:]=self.actmap
            if self.rotation==None or self.rotation==0:
                pass
            elif type(self.rotation) in [type(0),type(0.)]:
                self.mirrorSurface.rotate(self.rotation)
            elif type(self.rotation)!=type(None):
                self.mirrorSurface.rotate(self.rotation())
        else:
            #self.mirrorSurface.phsOut[:]=0#no need, since it doesn't get used.
            pass

    def applySlaving(self,actmap,actslaves):
        """actslaves is a dict of indx:slavelist) where indx is the index of the actuator to be slaved, and slavelist is a list of (indx,val) where indx is the actuator index to get slaving from, and val is the strength to apply.
        This allows any given actuator to be slaved to any combination of others.  It can also be slaved to itself...
        """
        for indx in actslaves.keys():
            slavelist=actslaves[indx]
            act=0.
            for pindx,val in slavelist:
                act+=actmap[pindx]*val
            actmap[indx]=act
        return actmap

    def getSlaves(self):
        """primarily for testing/gui use - plots all the slaved actuators
        """
        actmap=numpy.zeros((self.nact*self.nact,),numpy.int32)
        for indx in self.actSlaves.keys():
            actmap[indx]=1
        actmap.shape=self.nact,self.nact
        return actmap

    def makeMirrorModes(self,modeMatrix=None,fitpup=1):
        """Called by the GUI, after it has set up the modeMatrix object
        (which it will get acter a call to the makeMirrorModes object
        in the reconstructor)
        If called with modeMatrix==None, will do individual pokes.
        """
        actmx=numpy.zeros((self.nact,self.nact),numpy.float32)
        actmx1d=numpy.reshape(actmx,(self.nact*self.nact,))
        self.mirrorModes=numpy.zeros((self.nacts,self.npup,self.npup),numpy.float32)
        self.mirrorScale=numpy.zeros((self.nacts),numpy.float32)
        if self.pupil==None:
            raise Exception("Must specify a pupil function for making mirror modes")
        pfn=self.pupil.fn
        for i in range(self.nacts):
            if modeMatrix==None:
                y=self.dmindices[i]/self.nact
                x=self.dmindices[i]%self.nact
                actmx[y,x]=1.
            else:
                numpy.put(actmx1d,self.dmindices,modeMatrix[:,-1-i])
            tmp=self.mirrorSurface.fit(actmx)
            self.mirrorModes[i]=tmp
            if fitpup:
                self.mirrorModes[i]*=pfn
            self.mirrorScale[i]=numpy.sqrt(numpy.sum(self.mirrorModes[i]*self.mirrorModes[i]))
            self.mirrorModes[i]/=self.mirrorScale[i]

            if modeMatrix==None:
                actmx[y,x]=0.
        actmx[:]=0
        self.mirrorSurface.fit(actmx)
        util.FITS.Write(self.mirrorModes,"mirrorModes.fits")
        return self.mirrorModes

    def montePhaseCovariance(self):
        a=numpy.zeros((self.nacts,),numpy.float64)
        phs=self.parent["atmos"].outputData
        if type(self.mirrorModes)==type(None):
            self.makeMirrorModes(fitpup=self.fitPupMirrorModes)
        start=0
        for i in xrange(start,self.nacts):
            a[i]=numpy.sum(phs*self.mirrorModes[i])
        if type(self.monteNoll)==type(None):
            self.monteNoll=numpy.zeros((self.nacts,self.nacts),numpy.float64)
            self.navNoll=0
        self.monteNoll+=a[:,numpy.newaxis]*a[numpy.newaxis,:]
        self.navNoll+=1
    def finishMontePhaseCovariance(self,fname="montePhaseCovMat.fits"):
        self.montePhaseCovMat=self.monteNoll/(self.navNoll)#*self.pupil.sum)#*self.pxlscale**2)
                    # should pupil.sum be with or without the central obscuration? (not that
                    # it makes a huge difference!).  Don't scale to pixels, because the poke
                    # matrix scales from radians to pixels.
        if fname!=None:
            print "%s has been written"%fname
            util.FITS.Write(self.montePhaseCovMat,fname)
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
        if self.sentPlotsCnt==len(self.thisObjList)-1:
            txt+="""<plot title="xinterp_dm output%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        txt+="""<plot title="xinterp_dm mirror%s" cmd="data=-%s.dmphs*%s.thisObjList[%d].lineOfSight.wavelengthAdjustor" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,objname,self.sentPlotsCnt)
        txt+="""<plot title="xinterp_dm selected mirror%s" cmd="data=-%s.thisObjList[%d].lineOfSight.selectedDmPhs" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt)
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
        paramList.append(base.dataType.dataType(description="nAct",typ="eval",val="wfs_nsubx+1",comment="TODO: Number of actuators across the pupil"))
        paramList.append(base.dataType.dataType(description="wfs_n",typ="eval",val="npup/wfs_nsubx",comment="TODO: Number of pixels across a subap"))
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))
        paramList.append(base.dataType.dataType(description="dmObj",typ="code",val="import util.dm;dmObj=util.dm.dmOverview(dmInfoList,atmosGeom=atmosGeom)",comment="TODO: dmInfoList is a list of util.dm.dmInfo objects, initialised with (label (for this dm), idlist (list of (dm ID,source ID) or just sourceID, the idstr for a particular DM object (height and direction) and the idstr for a given source direction), height, nact, and various other things for which the defaults may be okay (see the source code))"))
        
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(npup,npup/2,npup/2*telSec/telDiam,wfs_nsubx,wfs_minarea)",comment="Telescope pupil object"))
        #paramList.append(base.dataType.dataType(description="sourceLam",typ="f",val="1650.",comment="source wavelength"))
        #paramList.append(base.dataType.dataType(description="wfslam",typ="f",val="1650.",comment="source wavelength"))
        #paramList.append(base.dataType.dataType(description="dmInterpType",typ="s",val="spline",comment="interpolation for DM surface"))

        return paramList
	

### Startup ############################################################
def makeMirrorMatrix(fout=None):
    """Make mirror matrix - ie the mirror shape for each poked actuator."""
    import sys
    sys.argv.append("--nostdin")
    class recon:
        outputData=None
        dataValid=1
    import util.Ctrl,util.FITS
    ctrl=util.Ctrl.Ctrl(globals=globals())
    paramfile=ctrl.paramfile
    idstr=ctrl.simID
    if idstr=="":
        idstr=None
    r=recon()
    this=dm(r,ctrl.config,args={},idstr=idstr)
    this.finalInitialisation()
    pokeval=ctrl.config.getVal("pokeval",default=1.)
    r.outputData=numpy.zeros((this.nacts,),numpy.float32)
    result=numpy.empty((this.nacts,this.npup*this.npup),numpy.float32)
    for i in xrange(this.nacts):
        r.outputData[i]=pokeval
        this.doNextIter()
        r.outputData[i]=0.
        result[i]=this.outputData.ravel()
    result/=pokeval
    if fout!=None:
        util.FITS.Write(result,fout)
    ctrl.go=0
    ctrl.mainloop([])
    return result

if __name__=="__main__":
    makeMirrorMatrix("tmp.fits")

