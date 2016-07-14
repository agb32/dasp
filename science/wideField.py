"""Module wideField - takes one or more iscrn outputs and combines together,
does interpolation etc along slightly different lines of sight, and from these, creates a wide field of view image, or images (sub-apertures).


Works as follows:

Generate a pupil phase.
Create the SHS images.
Convolve each with the correct part of the image.
Multiply by an enveloping cosine function.
Add into the final wide-field array.

Repeat for the next (slightly different) line of sight.


"""

import base.aobase,util.getNewCols
import numpy
import time,types
import iscrn
import wfscent
import util.atmos
import scipy.signal
#import util.dist,util.zernikeMod



class WideField(base.aobase.aobase):
    """Create a WideField object.  This object can take several iscrn
    objects as parents, and returns a wide-field image, or sub-aperture images of the relevant phase.

    Thoughts on resource sharing: each iatmos object stores phasescreens from
    the nLayers iscrn.thisObjDict objects.  This is potentially a large amount of data
    when replicated over a number of different iatmos objects (different sky
    directions).  However, the phasescreens should be identical.  So, we can
    create 1 iatmos object which resource shares between all the others.

    Alternatively, we have 1 iatmos object which is able to compute
    lots of directions from one call to generateNext.  You would then only
    ever need one iatmos module running on a given node, serving up pupil
    phase for lots of different directions, ie the outputData would be 3D.

    Having a resource sharing object is best - particularly if you have lots
    of directions - in the non-sharing case, outputData size is proportional to
    the number of direcions, while in the sharing case it isn't.

    Resource sharing implemented.


    Parameters of relevance:
    nFieldX, nFieldY - the number of directions evaluated when computing.  Typically, nfieldY isn't specified and so assumed same as nfieldX

    widefieldImage:  Array of size nsubx,nsubx,fftsize*(nfieldy+1.)/2,fftsize*(nfieldx+1.)/2
    The pixelscale of this image will be that of the unbinned image, i.e based on n and nfft.  ie util.calcPxlScale.pxlScale(lam,telDiam/nsubx,n,nfft,1)

    widefieldImageBoundary: Typically a small integer (or 0).  Try increasing this if there are obvious rectangles within the images - this increases the individual psf sizes.

    fov: Within the wfsOverview, a fov should be specified for the wfs.  This is equal to the lgspsf pixel size * lgs_nfft*(nfieldX-1.)/2.  *2.

    clipsize for this wfs should be equal to nfft.

    wideFieldDmList: A list of DM parents.  Optional - if not specified, these will be assumed to be parents starting with "dm".

    fieldAlt: None, or an array of shape nfieldY,nfieldX, containing the source heights for this part of the field (e.g. for differential cone effect using LGSs).

    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Initialise the object.  The parent object here should be a dictionary with values for each atmospheric layer, ie instances of iscrn, or DMs which are sending their full surface.  """
        if type(parent)!=type({}):
            parent={"L0":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.scrnDataType="d"#self.config.getVal("dataType")
        self.outDataType="f"
        if forGUISetup==1:
            obj=self.config.getVal("wfsOverview")
            self.wfsobj=obj.getWfsByID(idstr)
            self.nsubx=self.wfsobj.nsubx
            self.nimg=self.wfsobj.nimg
            self.nFieldX=self.config.getVal("nFieldX")
            self.nFieldY=self.config.getVal("nFieldY",default=self.nFieldX)
            self.outputData=[(self.nsubx,self.nsubx,(self.nFieldY-1)*self.nimg/2+self.nimg,(self.nFieldX-1)*self.nimg/2+self.nimg),self.outDataType]
        else:#setup for real
            self.parentDict=parent
            self.dmParentDict={}
            self.scrnParentDict={}
            dmParentList=self.config.getVal("wideFieldDmList",default=[])
            self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
            for key in self.parentDict.keys():
                try:
                    self.dmObj.getDM(key,raiseerror=1)
                    isdm=1
                except:#either not a DM, or no dmObj specified...
                    isdm=0
                if isdm or key.startswith("dm") or key in dmParentList:
                    self.dmParentDict[key]=self.parentDict[key]
                else:
                    self.scrnParentDict[key]=self.parentDict[key]
            print "widefield screens: %s"%str(self.scrnParentDict.keys())
            print "widefield dms: %s"%str(self.dmParentDict.keys())
            self.ndm=len(self.dmParentDict)
            self.dmKeys=self.dmParentDict.keys()#fix the order of parents -arbitrary
            self.npup=self.config.getVal("npup")#pupil size (pixels)
            self.timing=self.config.getVal("timing",default=0)
            self.interpolationNthreads=self.config.getVal("interpolationNthreads",default=0)#a tuple of (nthreads,nblockx,nblocky)
            self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
            self.dmList=[]
            for key in self.dmParentDict.keys():
                self.dmList.append(self.dmObj.getDM(key,raiseerror=1))

            self.parentSendWholeScreen=0
            self.doneFinalInit=0
            self.sentPlotsCnt=0
            self.layerListDict={}
            self.pupilphs=numpy.zeros((self.npup,self.npup),numpy.float32)
            self.atmosGeom=self.config.getVal("atmosGeom")
            self.nFieldX=self.config.getVal("nFieldX")
            self.nFieldY=self.config.getVal("nFieldY",default=self.nFieldX)
            self.widefieldImage=self.config.getVal("widefieldImage")
            #If there are visible boundaries in the image then increase the widefieldImageBoundary (and adjust the widefieldImage to suit).
            self.boundary=self.config.getVal("widefieldImageBoundary",default=0)
            self.smoothingPsf=self.config.getVal("smoothingPsf",raiseerror=0)
            self.fieldAlt=self.config.getVal("fieldAlt",raiseerror=0)
            if self.fieldAlt is not None:
                if self.fieldAlt.shape!=(self.nFieldY,self.nFieldX):
                    raise Exception("fieldAlt must be of shape %d,%d"%(self.nFieldY,self.nFieldX))
            obj=self.config.getVal("wfsOverview")
            self.wfsobj=obj.getWfsByID(idstr)
            self.nsubx=self.wfsobj.nsubx
            self.nimg=self.wfsobj.nimg
            self.fftsize=self.wfsobj.nfft
            self.clipsize=self.wfsobj.clipsize#This determines the shift... 
            if self.clipsize!=self.fftsize:
                raise Exception("wideField wfs objects must not be clipped - i.e. clipsize should equal fftsize")#This is because we have to add the PSFs in their relative positions.
            #We add the noise at the end, not within the wfscent objects...
            self.addPoisson=self.wfsobj.addPoisson
            self.readoutNoise=self.wfsobj.readoutNoise
            self.bglevel=self.wfsobj.bglevel
            self.wfsobj.addPoisson=0
            self.wfsobj.readoutNoise=0.
            self.wfsobj.bglevel=0.

            self.psfsize=self.fftsize*2#to avoid aliasing
            self.cos2d=numpy.zeros((self.nimg,self.nimg),numpy.float32)#not actually cos... quadratic.
            lin=1-numpy.abs(numpy.arange(self.nimg)-(self.nimg-1)/2.)/((self.nimg-1)/2.)
            self.cos2d[:]=lin
            self.cos2d=self.cos2d*self.cos2d.T#When added together half a spacing apart, this function is flat.
            #And now we need these window functions for corner and edge cases.
            self.cos2dBottom=self.cos2d.copy()
            self.cos2dBottom[:self.nimg/2]+=self.cos2d[self.nimg/2:]
            self.cos2dTop=self.cos2d.copy()
            self.cos2dTop[-self.nimg/2:]+=self.cos2d[:self.nimg/2]
            self.cos2dLeft=self.cos2d.copy()
            self.cos2dLeft[:,:self.nimg/2]+=self.cos2d[:,self.nimg/2:]
            self.cos2dRight=self.cos2d.copy()
            self.cos2dRight[:,-self.nimg/2:]+=self.cos2d[:,:self.nimg/2]
            
            self.cos2dBottomLeft=self.cos2dBottom.copy()
            self.cos2dBottomLeft[:,:self.nimg/2]+=self.cos2dBottom[:,-self.nimg/2:]
            self.cos2dBottomRight=self.cos2dBottom.copy()
            self.cos2dBottomRight[:,-self.nimg/2:]+=self.cos2dBottom[:,:self.nimg/2]
            self.cos2dTopLeft=self.cos2dTop.copy()
            self.cos2dTopLeft[:,:self.nimg/2]+=self.cos2dTop[:,-self.nimg/2:]
            self.cos2dTopRight=self.cos2dTop.copy()
            self.cos2dTopRight[:,-self.nimg/2:]+=self.cos2dTop[:,:self.nimg/2]

            #Also cases where nFieldX or Y is 1 - eg maybe for an elongated sodium spot...
            self.cos2dTopBottom=self.cos2dBottom.copy()
            self.cos2dTopBottom[-self.nimg/2:]+=self.cos2d[:self.nimg/2]
            self.cos2dLeftRight=self.cos2dLeft.copy()
            self.cos2dLeftRight[:,-self.nimg/2:]+=self.cos2d[:,:self.nimg/2]
            self.cos2dTopBottomLeft=self.cos2dTopBottom.copy()
            self.cos2dTopBottomLeft[:,:self.nimg/2]+=self.cos2dTopBottom[:,self.nimg/2:]
            self.cos2dTopBottomRight=self.cos2dTopBottom.copy()
            self.cos2dTopBottomRight[:,-self.nimg/2:]+=self.cos2dTopBottom[:,:self.nimg/2]
            self.cos2dLeftRightTop=self.cos2dLeftRight.copy()
            self.cos2dLeftRightTop[-self.nimg/2:]+=self.cos2dLeftRight[:self.nimg/2]
            self.cos2dLeftRightBottom=self.cos2dLeftRight.copy()
            self.cos2dLeftRightBottom[:self.nimg/2]+=self.cos2dLeftRight[self.nimg/2:]

            self.studySubap=self.config.getVal("studySubap",(0,0))#the subap of interest, for debugging purposes.  Can be changed from the GUI to view a different one.
            self.studyPsfs=numpy.zeros((self.nFieldY,self.nFieldX,self.fftsize+2*self.boundary,self.fftsize+2*self.boundary),numpy.float32)
            self.studySHS=numpy.zeros((self.nFieldY,self.nFieldX,self.nimg,self.nimg),numpy.float32)
            
            self.psf=numpy.zeros((self.nsubx,self.nsubx,self.psfsize,self.psfsize),numpy.float32)
            self.interpNthreads=self.config.getVal("interpolationNthreads",default=0
)

            self.directPhaseScreen=self.config.getVal("directPhaseScreen",default=1)#are we allowed to access parent.screen if parent is an iscrn object?  
            self.telDiam=self.config.getVal("telDiam")
            self.ntel=self.config.getVal("ntel",default=self.npup)#tel diameter in pxls.  
            self.pupil=self.config.getVal("pupil")
            #self.telSec=self.config.getVal("telSec")
            self.scrnScale=self.telDiam/float(self.ntel)#self.config.getVal("scrnScale")
            scrnScale=self.scrnScale
            self.pixScale=self.telDiam/float(self.npup)
            self.jmax=self.config.getVal("jmax",default=55)#jmax here is by default 55.  A better value might be npup*0.75 ^2/2, about the most you would hope.  If npup is large, this will take a long time, so you should probably use less for jmax.  Only used for zernike variance (from GUI).
            self.zernObj=None

            self.windDirection=self.atmosGeom.windDirections()
            self.vWind=self.atmosGeom.windSpeeds()
            self.layerAltitude=self.atmosGeom.layerAltitudes()#these are scaled by zenith.
            self.r0=self.atmosGeom.r0
            self.L0=self.atmosGeom.l0
    
            self.control={"cal_source":0,"fullPupil":0,"removePiston":0}#full pupil is used as a flag for phase profiling (xinterp_recon etc) if want to return the whole pupil.  The removePiston is used somewhere too.  
            self.tstep=self.config.getVal("tstep")
            #self.niters=0
            self.outputData=numpy.zeros((self.nsubx,self.nsubx,(self.nFieldY-1)*self.nimg/2+self.nimg,(self.nFieldX-1)*self.nimg/2+self.nimg),self.outDataType)#resource sharing
            #self.scrnXPxls=self.config.getVal("scrnXPxls")#will depend on which screen it is, so this is no good!
            #self.scrnYPxls=self.config.getVal("scrnYPxls")
            #Get the initial phase screen.
            self.phaseScreens={}

            # The initial phase (not passed from iscrn) should be created and passed into this object.  Copy and convert to numpy array (from numarray)
            #Note, if you are allowing directPhaseScreens, and wish not to bother creating an initial phase, use a dummy object in initPhs, with a copy method and a shape tuple.  This should then work saving you time and memory.
            self.phStructFn=None
            self.phStructFnIter=0
            self.coeffsZernike=None
            self.coeffsZernike2=None
            self.zernikeIters=0
            if self.config.getVal("precomputeYGradient",default=0):#set to 1 for faster operation with larger memory comsumption
                self.ygradient={}
            else:
                self.ygradient=None
            for pkey in self.scrnParentDict.keys():
                #each iscrn parent can have 1 or many layers.  So, get this into a dictionary.
                self.layerListDict[pkey]=self.config.getVal("layerList",{},searchOrder=["iscrn_%s","iscrn","globals"]).get(pkey,[pkey])


                for key in self.layerListDict[pkey]:
                    #First, see if we can copy from parent - if not, create ourselves (this relies on the seed being a constant, not the current time).
                    try:
                        self.phaseScreens[key]=self.scrnParentDict[pkey].thisObjDict[key].screen#[xxx].copy()
                        print "iatmos: Copied initial screen from parent %s[%d]"%(str(pkey),str(key))
                        print "TODO: iatmos - is this copy of initial screen needed?"
                    except:
                        print "iatmos: cannot copy parent screen %s - generating directly."%str(key)
                        self.phaseScreens[key]=numpy.array(iscrn.computeInitialScreen(self.config,idstr=str(key)))
                    if self.phaseScreens[key].dtype.char!=self.scrnDataType:
                        raise Exception("iatmos and iscrn should use the same dataType value %s %s"%(self.phaseScreens[key].dtype.char,self.scrnDataType))
                    ps=self.phaseScreens[key]
                    if ps.shape!=(self.atmosGeom.getScrnYPxls(key,rotateDirections=1),self.atmosGeom.getScrnXPxls(key,rotateDirections=1)):
                        raise Exception("Phase screen size unexpected in iatmos: %s -> %s"%(str(ps.shape),str((self.atmosGeom.getScrnYPxls(key,rotateDirections=1),self.atmosGeom.getScrnXPxls(key,rotateDirections=1)))))
                    if self.ygradient!=None:
                        self.ygradient[key]=numpy.empty(ps.shape,self.scrnDataType)
                        self.ygradient[key][1:-1]=(ps[2:]-ps[:-2])*0.5
                        self.ygradient[key][0]=ps[1]-ps[0]
                        self.ygradient[key][-1]=ps[-1]-ps[-2]

            #Now get the layer offsets which are used to determine where in a layer to start looking for this source.  Basically, you can use sourceTheta and sourcePhi and sourceAlt to determine the distance in m from the on axis location.  Convert this into pixels, and you may get a negative number (ie if it is to the left of the onaxis position) depending on which source and which layer.  The maximum possible negative number is then equal to this layerOffset, which makes sure that all indexes into the array are positive. 
            self.layerXOffset=self.atmosGeom.getLayerXOffset(rotateDirections=1)
            self.layerYOffset=self.atmosGeom.getLayerYOffset(rotateDirections=1)
            self.scrnXPxls={}
            self.scrnYPxls={}
            #self.interpPhs=numpy.zeros((self.npup,self.npup),self.dataType)#resource sharing
            for pkey in self.scrnParentDict.keys():
                for key in self.layerListDict[pkey]:
                    self.scrnXPxls[key]=self.phaseScreens[key].shape[1]#agbc from 0
                    self.scrnYPxls[key]=self.phaseScreens[key].shape[0]#agbc from 1
            self.rowAdd={}
            self.maxRowAdd={}
            self.newRows={}
            self.insertPos={}
            #create objects for working out the interpolation required for each new iteration...
            for pkey in self.scrnParentDict.keys():
                for key in self.layerListDict[pkey]:
                    self.rowAdd[key]=self.vWind[key]/self.pixScale*self.tstep#number of pixels to step each iteration (as float).  
                    self.newRows[key]=util.getNewCols.getNewCols(numpy.fabs(self.rowAdd[key]))
                    self.maxRowAdd[key]=int(numpy.ceil(numpy.fabs(self.rowAdd[key])))
                    self.insertPos[key]=0
            # we add new rows at the bottom of the array (note, that if plotting in Gist, this is the top of the array).
            self.inputData={}#hold the new phase cols/rows.
            self.rowInput={}
            self.interpPosRow={}
            self.nremRow={}
            self.naddRow={}
            for key in self.scrnParentDict.keys():
                self.inputData[key]=None
                self.rowInput[key]=None
            self.initialise(self.parent,self.idstr[0])
            #for i in xrange(len(self.idstr)):
            #    idstr=self.idstr[i]
            #    parent=self.parentList[i]
            #    self.initialise(parent,idstr)
    def newParent(self,parent,idstr=None):
        raise Exception("iatmos - not yet able to accept new parent... (needs some extra coding)")

    def initialise(self,parentDict,idstr):
        """note, parent should be a dictionary here...
        This is called both from the __init__ method, and from the addNewIdObject method.  

        """
        if type(parentDict)!=type({}):#expecting a dictionary of parents
            parentDict={"L0":parentDict}
            print "WARNING - wideField: I was expecting a dictionary for parent here..."
        if not set(parentDict)==set(self.parentDict):#check they have teh same parents (main and resource sharer).
            raise Exception("wideField - all resource sharing objects must have same parents... (though whether they use them can be selected in the param file)")
        this=base.aobase.resourceSharer(parentDict,self.config,idstr,self.moduleName)
        self.thisObjList.append(this)
        
        sourceTheta=self.atmosGeom.sourceTheta(idstr)
        sourceAlt=self.atmosGeom.sourceAlt(idstr)
        sourcePhi=self.atmosGeom.sourcePhi(idstr)
        sourceLam=self.atmosGeom.sourceLambda(idstr)
        zenith=self.atmosGeom.zenith
        fov=self.atmosGeom.sourceFov(idstr)
        if sourceLam==None:
            sourceLam=self.config.getVal("sourceLam")
            print "Warning - atmosGeom sourceLam==None, using %g"%sourceLam
        layerList=this.config.getVal("atmosLayerList",default=None,raiseerror=0)
        if layerList==None:
            layerList=[]
            for key in self.layerListDict.keys():
                layerList+=self.layerListDict[key]
        intrinsicPhase=this.config.getVal("intrinsicPhase",default=None,raiseerror=0)
        storePupilLayers=this.config.getVal("storePupilLayers",default=0)
        computeUplinkTT=this.config.getVal("computeUplinkTT",default=0)
        launchApDiam=this.config.getVal("launchApDiam",default=0)

        this.atmosObj=util.atmos.iatmos(sourceAlt,sourceLam,sourceTheta,sourcePhi,self.npup,self.pupil,self.rowAdd,self.layerAltitude,self.windDirection,self.phaseScreens,self.ygradient,self.scrnScale,self.layerXOffset,self.layerYOffset,layerList,zenith,intrinsicPhase=intrinsicPhase,storePupilLayers=storePupilLayers,computeUplinkTT=computeUplinkTT,launchApDiam=launchApDiam,ntel=self.ntel,telDiam=self.telDiam,interpolationNthreads=self.interpolationNthreads,outputData=self.pupilphs,fov=fov)

        class Parent:
            dataValid=1
            outputData=self.pupilphs
        self.wfsParent=Parent()
        this.wfsObj=wfscent.wfscent(self.wfsParent,self.config,idstr=idstr)
        this.wfsObj.finalInitialisation()

        sourceAlt=self.atmosGeom.sourceAlt(idstr)
        sourceLam=self.atmosGeom.sourceLambda(idstr)
        this.lineOfSight=[]
        #interpolate=numpy.zeros((self.npup,self.npup),numpy.float32)
        pitchY=this.atmosObj.fov*2*2/(self.nFieldY+1.)#fov sampling pitch.
        pitchX=this.atmosObj.fov*2*2/(self.nFieldX+1.)
        sourceTheta=self.atmosGeom.sourceTheta(idstr)#*numpy.pi/180./3600.
        sourcePhi=self.atmosGeom.sourcePhi(idstr)#*numpy.pi/180.
        for n in range(self.ndm):
            dm=self.dmList[n]
            key=self.dmKeys[n]
            wavelengthAdjustor=dm.reconLam/float(sourceLam)
            dmpup=dm.calcdmpup(self.atmosGeom)
            for i in range(self.nFieldY):
                for j in range(self.nFieldX):
                    ydiff=-this.atmosObj.fov+pitchY/2.+i*pitchY/2.
                    xdiff=-this.atmosObj.fov+pitchX/2.+j*pitchX/2.

                    xcentre=sourceTheta*numpy.cos(sourcePhi*numpy.pi/180.)
                    ycentre=sourceTheta*numpy.sin(sourcePhi*numpy.pi/180.)
                    xnew=xcentre+xdiff
                    ynew=ycentre+ydiff
                    thetaNew=numpy.sqrt(xnew*xnew+ynew*ynew)*numpy.pi/180./3600.
                    phiNew=numpy.arctan2(ynew,xnew)
                    los=util.dm.DMLineOfSight(dmpup,self.npup,dm.height,self.dmParentDict[key].outputData,sourceAlt,thetaNew,phiNew,self.telDiam,wavelengthAdjustor,nthreads=self.interpNthreads)#,interpolated=interpolate)
                    this.lineOfSight.append(los)


    def finalInitialisation(self):
        """since all memories are the same size (npup can't change...), its
        easy here..."""
        if self.doneFinalInit:
            return
        self.doneFinalInit=1
        for this in self.thisObjList:
            this.atmosObj.initMem()#self.outputData)#,self.interpPhs)

    def generateNext(self,msg=None):
        """
        This function is called when it is okay to produce the next iteration
        of the simulation.
        Not expecting any msgs.
        """
        t1=time.time()
        atmosObj=self.thisObjList[self.currentIdObjCnt].atmosObj
        sourceTheta=atmosObj.sourceTheta
        sourcePhi=atmosObj.sourcePhi
        wfsObj=self.thisObjList[self.currentIdObjCnt].wfsObj
        if self.debug:
            print "iatmos: GenerateNext (debug=%s)"%str(self.debug)
        if self.generate==1:
            if self.newDataWaiting:
                nvalid=0
                if self.currentIdObjCnt==0:
                    #if its the first object, we assume all screens are ready
                    #and construct them.   
                    for pkey in self.scrnParentDict.keys():
                        if self.scrnParentDict[pkey].dataValid==1:
                            if self.inputData[pkey] is not self.scrnParentDict[pkey].outputData:
                                #create the arrays...
                                #print "Allocating iatmos arrays (hopefully this only happens during first iteration...)"
                                self.inputData[pkey]=self.scrnParentDict[pkey].outputData
                                if self.inputData[pkey].dtype.char!=self.scrnDataType:
                                    raise Exception("iatmos and iscrn dataTypes must be same (%s %s, key %s)"%(self.inputData[pkey].dtype.char,self.scrnDataType,pkey))
                                if self.parentSendWholeScreen==0:
                                    pos=0
                                    for key in self.layerListDict[pkey]:
                                        self.rowInput[key]=self.inputData[pkey][pos:pos+self.maxRowAdd[key]*self.scrnXPxls[key]]
                                        pos+=self.maxRowAdd[key]*self.scrnXPxls[key]
                                        self.rowInput[key].shape=(self.maxRowAdd[key],self.scrnXPxls[key])
                            nvalid+=1
                    if nvalid==len(self.scrnParentDict.keys()):
                        self.dataValid=1
                        self.makeLayers()#this is the first of resource sharers
                    elif nvalid>0:
                        print "ERROR: iatmos - received wrong number (%d/%d) of phase screens"%(self.nvalid,len(self.scrnParentDict.keys()))
                        self.dataValid=0
                    else:
                        print "iatmos: Waiting for data from iscrn, but not valid"
                        self.dataValid=0
                if self.dataValid:
                    self.outputData[:]=0
                    fovpitchY=atmosObj.fov*2*2/(self.nFieldY+1.)
                    fovpitchX=atmosObj.fov*2*2/(self.nFieldX+1.)
                    altOrig=atmosObj.sourceAlt
                    thetaOrig=atmosObj.sourceTheta
                    phiOrig=atmosObj.sourcePhi
                    for fieldY in range(self.nFieldY):
                        for fieldX in range(self.nFieldX):
                            if self.control["cal_source"]:#calibration source
                                self.pupilphs[:]=0.
                            else:
                                #update for field
                                #We have nFieldX/Y covering this field, and each sub-field overlaps by 50%.  Therefore, (nFieldX+1)/2*subfov - which should equal the full fov.  
                                #ydiff=-atmosObj.fov+pitchY/2.+fieldY*pitchY/2.
                                #xdiff=-atmosObj.fov+pitchX/2.+fieldX*pitchX/2.
                                ydiff=(-(self.nFieldY-1)/2.+fieldY)*fovpitchY/2.
                                xdiff=(-(self.nFieldX-1)/2.+fieldX)*fovpitchX/2.


                                #ydiff=(fieldY/float(self.nFieldY-1)-0.5)*atmosObj.fov
                                #xdiff=(fieldX/float(self.nFieldX-1)-0.5)*atmosObj.fov
                                xcentre=sourceTheta*numpy.cos(sourcePhi*numpy.pi/180.)
                                ycentre=sourceTheta*numpy.sin(sourcePhi*numpy.pi/180.)
                                xnew=xcentre+xdiff
                                ynew=ycentre+ydiff
                                thetanew=numpy.sqrt(xnew*xnew+ynew*ynew)
                                phinew=numpy.arctan2(ynew,xnew)
                                atmosObj.sourceTheta=thetanew
                                atmosObj.sourcePhi=phinew
                                if self.fieldAlt!=None:
                                    atmosObj.sourceAlt=self.fieldAlt[fieldY,fieldX]
                                #print "%g diff %g, %g, pitch %g, %g, dir: %g, %g"%(atmosObj.fov, xdiff,ydiff,fovpitchX,fovpitchY,thetanew,phinew)
                                atmosObj.updatePositionDict()
                                atmosObj.createPupilPhs(self.interpPosRow,self.insertPos,self.control)
                            #Add in any DM phase... (select correct part for part of the fov)
                            self.addDmPhase(self.pupilphs,fieldX,fieldY)
                            #Now, use the phase to generate the wide-field images...
                            #update the psf for correct part of the field.
                            #widefieldImageFov tells us which parts to select.  
                            ny=self.widefieldImage.shape[2]-2*self.boundary
                            nx=self.widefieldImage.shape[3]-2*self.boundary
                            if nx<self.fftsize*(self.nFieldX+1)/2:
                                raise Exception("Wrong size widefieldImage - needs to be at least %d."%(self.fftsize*(self.nFieldX+1)/2))
                            pitchX=nx//(self.nFieldX+1)#and is this many pixels - the +1 allows for edge effects.
                            pitchY=ny//(self.nFieldY+1)#and is this many pixels - the +1 allows for edge effects.
                            startX=int((nx-self.fftsize*(self.nFieldX+1.)/2)/2)
                            startY=int((ny-self.fftsize*(self.nFieldY+1.)/2)/2)
                            #print self.psf.shape,self.widefieldImage.shape,startY+fieldY*pitch,startY+fieldY*pitch+self.fftsize,startX+fieldX*pitch,startX+fieldX*pitch+self.fftsize
                            #print fieldY,fieldX,self.widefieldImage.shape,startY,pitchY,self.fftsize,startX,pitchX
                            for y in range(self.nsubx):
                                for x in range(self.nsubx):
                                    self.psf[y,x,self.fftsize/2-self.boundary:self.fftsize/2*3+self.boundary,self.fftsize/2-self.boundary:self.fftsize/2*3+self.boundary]=self.widefieldImage[y,x,startY+fieldY*pitchY:startY+fieldY*pitchY+self.fftsize+2*self.boundary,startX+fieldX*pitchX:startX+fieldX*pitchX+self.fftsize+2*self.boundary]
                            #Store the psfs for 1 subap.
                            self.studyPsfs[fieldY,fieldX]=self.widefieldImage[self.studySubap[0],self.studySubap[1],startY+fieldY*pitchY:startY+fieldY*pitchY+self.fftsize+2*self.boundary,startX+fieldX*pitchX:startX+fieldX*pitchX+self.fftsize+2*self.boundary]
                            wfsObj.thisObjList[0].wfscentObj.updatePsf(self.psf,updateSig=1)
                            #print numpy.any(numpy.isnan(self.pupilphs))
                            wfsObj.doNextIter()
                            wfsObj.currentIdObjCnt-=1#could get to below 0, but at the moment,that doesn't matter.
                            self.studySHS[fieldY,fieldX]=wfsObj.outputData[self.studySubap[0],self.studySubap[1]]
                            if wfsObj.dataValid:
                                #Now put the PSFs into the array.
                                self.prepareOutput(wfsObj.outputData,fieldX,fieldY)
                            else:#integration not yet ready.
                                raise Exception("This shouldn't happen - since we're using wfscent for many different directions")
                    atmosObj.sourceAlt=altOrig
                    atmosObj.sourceTheta=thetaOrig
                    atmosObj.sourcePhi=phiOrig

                    wfsObj.currentIdObjCnt+=1#ok
                    if self.smoothingPsf!=None:
                        for i in range(self.nsubx):
                            for j in range(self.nsubx):
                                self.outputData[i,j]=scipy.signal.fftconvolve(self.outputData[i,j],self.smoothingPsf,mode="same")
                    if self.control["cal_source"]==0:
                        if self.addPoisson!=0:
                            self.outputData[:]=numpy.random.poisson(self.outputData)
                        if self.readoutNoise!=0:
                            self.outputData+=numpy.random.normal(self.bglevel,self.readoutNoise)
                        elif self.bglevel!=0:
                            self.outputData+=self.bglevel

            else:#no new data ready
                self.dataValid=0
        else:
            self.dataValid=0
        #restore to central position
        atmosObj.sourceTheta=sourceTheta
        atmosObj.sourcePhi=sourcePhi
        atmosObj.updatePositionDict()#not sure we need to do this one, but anyway...

        if self.debug:
            print "wideField: done generateNext (debug=%s)"%str(self.debug)
        self.generateNextTime=time.time()-t1

    def addDmPhase(self,pupilphs,fieldX,fieldY):
        for i in range(self.ndm):
            #select the line of sight object for this DM.
            los=self.thisObjList[self.currentIdObjCnt].lineOfSight[i*self.nFieldX*self.nFieldY + fieldY*self.nFieldX + fieldX]
            
            #Add the correct portion of the DM to the pupilphs
            los.selectSubPupil(pupilphs,addToOutput=1,removePiston=0)

        
    def prepareOutput(self,wfsimg,fieldX,fieldY):
        ystep=self.cos2d.shape[0]/2
        xstep=self.cos2d.shape[1]/2
        sy=fieldY*ystep
        ey=fieldY*ystep+self.cos2d.shape[0]
        sx=fieldX*xstep
        ex=fieldX*xstep+self.cos2d.shape[1]
        fy=(wfsimg.shape[2]-self.cos2d.shape[0])//2
        fx=(wfsimg.shape[3]-self.cos2d.shape[1])//2
        if fieldX==0:
            if self.nFieldX==1:
                if fieldY==0:
                    if self.nFieldY==1:
                        window=0*self.cos2d+1#just 1... but why would you do a widefield with nFieldX/Y==1???
                        
                    else:
                        window=self.cos2dLeftRightBottom
                elif fieldY==self.nFieldY-1:
                    window=self.cos2dLeftRightTop
                else:
                    window=self.cos2dLeftRight
            else:
                if fieldY==0:
                    if self.nFieldY==1:
                        window=self.cos2dTopBottomLeft
                    else:
                        window=self.cos2dBottomLeft
                elif fieldY==self.nFieldY-1:
                    window=self.cos2dTopLeft
                else:
                    window=self.cos2dLeft
        elif fieldX==self.nFieldX-1:
            if fieldY==0:
                if self.nFieldY==1:
                    window=self.cos2dTopBottomRight
                else:
                    window=self.cos2dBottomRight
            elif fieldY==self.nFieldY-1:
                window=self.cos2dTopRight
            else:
                window=self.cos2dRight
        else:
            if fieldY==0:
                if self.nFieldY==1:
                    window=self.cos2dTopBottom
                else:
                    window=self.cos2dBottom
            elif fieldY==self.nFieldY-1:
                window=self.cos2dTop
            else:
                window=self.cos2d
        for y in range(self.nsubx):
            for x in range(self.nsubx):
                img=wfsimg[y,x]
                if img.shape!=window.shape:
                    #select the central part
                    img=img[fy:fy+window.shape[0],fx:fx+window.shape[1]]

                self.outputData[y,x,sy:ey,sx:ex]+=img*window


    def makeImage(self):
        """For simctrl.py use - i.e. for GUI"""
        od=self.outputData
        img=numpy.zeros((od.shape[0]*od.shape[2],od.shape[1]*od.shape[3]),numpy.float32)
        for i in range(od.shape[0]):
            for j in range(od.shape[1]):
                img[i*od.shape[2]:(i+1)*od.shape[2],j*od.shape[3]:(j+1)*od.shape[3]]=od[i,j]
        return img

    def makeLayers(self):
        for pkey in self.scrnParentDict.keys():#for each atmosphere layer...
            if self.directPhaseScreen and hasattr(self.scrnParentDict[pkey],"thisObjDict"):
                #we can just use the screens from the parents (though we'll still have to create the gradients).
                for key in self.layerListDict[pkey]:
                    
                    #use the exact same copy from the iscrn object.
                    #Note, this only works if it is in the same python process.
                    self.phaseScreens[key]=self.scrnParentDict[pkey].thisObjDict[key].screen
                    nremRow,naddRow,self.interpPosRow[key]=self.newRows[key].next()
                    self.insertPos[key]+=naddRow#work out the new insert position
                    if self.insertPos[key]>=self.scrnYPxls[key]:#wrap.
                        self.insertPos[key]-=self.scrnYPxls[key]
                    ps=self.phaseScreens[key]
                    #and update the gradients
                    for i in range(self.maxRowAdd[key]):
                        pos=self.insertPos[key]-self.maxRowAdd[key]+i
                        gpos=pos-1#gradient position
                        ggpos=pos-2
                        if pos<0:
                            pos+=self.scrnYPxls[key]#wrap
                        if gpos<0:
                            gpos+=self.scrnYPxls[key]#wrap
                        if ggpos<0:
                            ggpos+=self.scrnYPxls[key]#wrap
                        #ps[pos]=self.rowInput[key][i]
                        #and the gradient
                        if self.ygradient!=None:
                            self.ygradient[key][gpos]=(ps[pos]-ps[ggpos])*0.5
                    #and the last (partial) one...
                    if self.ygradient!=None:
                        self.ygradient[key][pos]=ps[pos]-ps[gpos]
            else:#construct the layer.
                for key in self.layerListDict[pkey]:
                    self.makeLayer(key)


    def makeLayer(self,key):
        """make the phasescreen from the newly added rows passed from iscrn object"""
        #first, add the new phase values into the phase array.
        nremRow,naddRow,self.interpPosRow[key]=self.newRows[key].next()
        self.nremRow[key]=nremRow
        self.naddRow[key]=naddRow
        if self.parentSendWholeScreen:
            #Note, this may raise an exception, and is probably an edge case anyway.
            self.phaseScreens[key]=self.scrnParentDict[key].outputData
        else:
            #here, we add newly created phase from iscrn to the existing phasescreen stored here.
            self.insertPos[key]+=naddRow#work out the new insert position
            if self.insertPos[key]>=self.scrnYPxls[key]:#wrap.
                self.insertPos[key]-=self.scrnYPxls[key]
            ps=self.phaseScreens[key]
            #and now insert up to that position.
            for i in range(self.maxRowAdd[key]):
                pos=self.insertPos[key]-self.maxRowAdd[key]+i
                gpos=pos-1#gradient position
                ggpos=pos-2
                if pos<0:
                    pos+=self.scrnYPxls[key]#wrap
                if gpos<0:
                    gpos+=self.scrnYPxls[key]#wrap
                if ggpos<0:
                    ggpos+=self.scrnYPxls[key]#wrap
                ps[pos]=self.rowInput[key][i]
                #and the gradient
                if self.ygradient!=None:
                    self.ygradient[key][gpos]=(ps[pos]-ps[ggpos])*0.5
            #and the last (partial) one...
            if self.ygradient!=None:
                self.ygradient[key][pos]=ps[pos]-ps[gpos]

    def tile(self,img):
        s=img.shape
        out=numpy.zeros((s[0]*s[2],s[1]*s[3]),numpy.float32)
        for i in range(s[0]):
            for j in range(s[1]):
                out[i*s[2]:i*s[2]+s[2],j*s[3]:j*s[3]+s[3]]=img[i,j]
        return out

    
    def createSingleLayerPhase(self,key):
        """Primarily for GUI... - gets the phase for this pupil at one height."""
        data=self.thisObjList[self.currentIdObjCnt].atmosObj.createSingleLayerPhs(self.interpPosRow,self.insertPos,key,self.control)
        return data

    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=iatmos.iatmos(...) then $OBJ would be replaced by scrn."""
        this=self.thisObjList[self.sentPlotsCnt]
        if this.idstr==None or this.idstr=="":
            id=""
        else:
            id=" (%s)"%this.idstr
        txt="""<plot title="output data" cmd="data=%s.makeImage()" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(objname)
        txt+="""<plot title="PSFs of 1 subap" cmd="data=%s.studyPsfs" ret="data" type="pylab" when="cmd" palette="gray"/>"""%(objname)
        txt+="""<plot title="Tiled PSFs of 1 subap" cmd="data=%s.tile(%s.studyPsfs)" ret="data" type="pylab" when="cmd" palette="gray"/>"""%(objname,objname)
        txt+="""<plot title="SHS Images of 1 subap" cmd="data=%s.studySHS" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(objname)
        txt+="""<plot title="Tiled SHS Images of 1 subap" cmd="data=%s.tile(%s.studySHS)" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(objname,objname)
        txt+="""<plot title="psf (of single direction)" cmd="data=%s.psf" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(objname)
        txt+="""<plot title="widefield psf" cmd="data=%s.widefieldImage" ret="data" type="pylab" when="cmd" palette="gray"/>"""%(objname)
        txt+="""<plot title="pupil phase (last position)" cmd="data=%s.pupilphs" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(objname)
        #txt+="""<plot title="Pupil phase screen%s (last resource sharer)" cmd="data=%s.pupilphs" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
        if len(self.thisObjList)==1:
            for pkey in self.scrnParentDict.keys():
                for key in self.layerListDict[pkey]:
                    txt+="""<plot title="layer phase %s%s" cmd="data=%s.createSingleLayerPhase('%s')" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(key,id,objname,key)
        txt+=self.thisObjList[0].wfsObj.plottable("%s.thisObjList[0].wfsObj"%objname)
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
        paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="dataType",typ="eval",val="this.globals.fpDataType",comment="Array numpy data type"))
        paramList.append(base.dataType.dataType(description="timing",typ="i",val="0",comment="Timing information"))
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="ntel",typ="eval",val="this.globals.npup",comment="Pixels for telescope"))
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(this.globals.npup,this.globals.ntel/2,this.globals.ntel/2*this.globals.telSec/this.globals.telDiam,this.globals.wfs_nsubx,this.globals.wfs_minarea)",comment="Telescope pupil"))
        paramList.append(base.dataType.dataType(description="tstep",typ="f",val="0.005",comment="TODO: timestep."))        
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by iscrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (iatmos), theta, phi, alt, nsubx or None)."))
        return paramList
