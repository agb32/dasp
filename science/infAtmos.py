#"""Module infAtmos - takes one or more infScrn outputs and combines together,
#does interpolation etc and creates the telescope pupil phase.
#"""


import base.aobase,util.getNewCols
import numpy
#import Numeric
import time,types
import infScrn
#from cmod.interp import mxinterp,linearshift
import util.atmos
import util.dist,util.zernikeMod
from scipy.special import gamma,kv
#if Numeric.__version__!="24.2":
    #raise Exception("Needs Numeric 24.2")
def calcLayerOffset(scrnSize,thetas,phis,altitude,npup,ntel,telDiam):
    """FUnction to compute the starting points for each layer (dictionary).  This is used in the parameter file."""
    print "ERROR: Please use util.atmos.calcLayerOffset"
    arcsecrad=2*numpy.pi/360./3600.
    degrad=2*numpy.pi/360.
    layerXOffset={}
    layerYOffset={}
    for altKey in altitude.keys():
        xpos=[]
        ypos=[]
        for sourceKey in thetas.keys():
            xpos.append(numpy.tan(thetas[sourceKey]*arcsecrad)*numpy.cos(phis[sourceKey]*degrad)*altitude[altKey]/telDiam*ntel-npup/2.)#scrnSize[altKey][0]/2.)
            ypos.append(numpy.tan(thetas[sourceKey]*arcsecrad)*numpy.sin(phis[sourceKey]*degrad)*altitude[altKey]/telDiam*ntel-npup/2.)#scrnSize[altKey][1]/2.)
        minx=min(xpos)
        miny=min(ypos)
        #print "minx,miny %g %g"%(minx,miny)
        if minx<0.:#less than zero...
            layerXOffset[altKey]=int(numpy.ceil(numpy.fabs(minx)))
        else:
            layerXOffset[altKey]=-int(numpy.ceil(minx))
        if miny<0.:#less than zero...
            layerYOffset[altKey]=int(numpy.ceil(numpy.fabs(miny)))
        else:
            layerYOffset[altKey]=-int(numpy.ceil(miny))
    layerOffset=(layerXOffset,layerYOffset)
    print layerOffset
    return layerOffset


class infAtmos(base.aobase.aobase):
    """Create an infAtmos object.  This object can take several infScrn
    objects as parents, and returns the pupil phase for a given source

    Thoughts on resource sharing: each infAtmos object stores phasescreens from
    the nLayers infScrn objects.  This is potentially a large amount of data
    when replicated over a number of different infAtmos objects (different sky
    directions).  However, the phasescreens should be identical.  So, we can
    create 1 infAtmos object which resource shares between all the others.

    Alternatively, we have 1 infAtmos object which is able to compute
    lots of directions from one call to generateNext.  You would then only
    ever need one infAtmos module running on a given node, serving up pupil
    phase for lots of different directions, ie the outputData would be 3D.

    Having a resource sharing object is best - particularly if you have lots
    of directions - in the non-sharing case, outputData size is proportional to
    the number of direcions, while in the sharing case it isn't.

    So, this object needs to implement resource sharing.  We therefore need to
    make sure that the input from infScrn objects is only used once for a given
    simultion iteration.  How to do this?  Look at the data to see if its
    changed? - ok but not elegant.  Or count the number of times called and
    update every time we're at the first resource sharer?  Better.  
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Initialise the object.  The parent object here should be a dictionary with values for each atmospheric layer, ie instances of infScrn.  The args dictionary should contain a key "initialPhase", the value of which is then a dictionary with keys the same as parent keys, and values equal to the initial phase screen.  This can be created using infScrn.computeInitialPhase()."""
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.dataType="d"#self.config.getVal("dataType")
        self.npup=self.config.getVal("npup")#pupil size (pixels)
        if forGUISetup==1:
            self.outputData=[(self.npup,self.npup),self.dataType]
        else:#setup for real
            self.timing=self.config.getVal("timing",default=0)
            self.doneFinalInit=0
            self.sentPlotsCnt=0
            # deg - rad conversion
            #self.degRad=self.config.getVal("degRad")
            # arcsec - radian conversion
            #self.arcsecRad=self.config.getVal("arcsecRad")

            #self.telFov=self.config.getVal("telFov")
            # Field of view (rad)
            #self.fov=self.telFov*self.arcsecRad  
            #fov=self.fov
            self.atmosGeom=self.config.getVal("atmosGeom",default=None,raiseerror=0)
            self.directPhaseScreen=self.config.getVal("directPhaseScreen",default=0)#are we allowed to access parent.screen if parent is an infScrn object?  
            self.parentSendWholeScreen=self.config.getVal("sendWholeScreen",default=0)#Is the infScrn object sending the whole screen or just the new bits?  This is similar to the directPhaseScreen but also works across MPI connections...
            self.telDiam=self.config.getVal("telDiam")
            self.ntel=self.config.getVal("ntel")#tel diameter in pxls.  
            self.pupil=self.config.getVal("pupil")
            #self.telSec=self.config.getVal("telSec")
            self.scrnScale=self.telDiam/float(self.ntel)#self.config.getVal("scrnScale")
            scrnScale=self.scrnScale
            self.pixScale=self.telDiam/float(self.npup)
            self.jmax=self.config.getVal("jmax",default=55)#jmax here is by default 55.  A better value might be npup*0.75 ^2/2, about the most you would hope.  If npup is large, this will take a long time, so you should probably use less for jmax.  Only used for zernike variance (from GUI).
            self.zernObj=None
            self.parentKeyList=None
            #The number of parent objects determines the number of layers.
            #self.nLayers=len(self.parentList[0].keys())
            if self.atmosGeom==None:#depreciated
                self.windDirection=self.config.getVal("windDirection")#depreciated Layer wind direction - a dictionary (one layer for each parent)
                self.vWind=self.config.getVal("vWind")#dictionary depreciated
                self.layerAltitude=self.config.getVal("altitude")#dictionary depreciated
                self.L0=self.config.getVal("l0")
                self.r0=self.config.getVal("r0")
            else:
                self.windDirection=self.atmosGeom.windDirections()
                self.vWind=self.atmosGeom.windSpeeds()
                self.layerAltitude=self.atmosGeom.layerAltitudes()#these are scaled by zenith.
                self.r0=self.atmosGeom.r0
                self.L0=self.atmosGeom.l0
    
            self.control={"cal_source":0,"profilePhase":0,"fullPupil":0,"removePiston":1}#full pupil is used as a flag for phase profiling (xinterp_recon etc) if want to return the whole pupil.
            self.tstep=self.config.getVal("tstep")
            #self.niters=0
            self.outputData=numpy.zeros((self.npup,self.npup),self.dataType)#resource sharing
            #self.scrnXPxls=self.config.getVal("scrnXPxls")#will depend on which screen it is, so this is no good!
            #self.scrnYPxls=self.config.getVal("scrnYPxls")
            #Get the initial phase screen.
            self.phaseScreens={}
            if not args.has_key("initialPhase"):
                print "Warning: infAtmos initial phases not found - using parent.keys(): %s"%str(self.parentList[0].keys())
                initPhs={}
                for key in self.parentList[0].keys():
                    initPhs[key]=""
            else:
                initPhs=args["initialPhase"]
                #print "WARNING: infAtmos initial phases not found.  This will mean the atmosphere is flat for the first few iterations.  You need to use infScrn.computeInitialScreen() method first."
                #raise Exception("Initial phases not found")
            # The initial phase (not passed from infScrn) should be created and passed into this object.  Copy and convert to numpy array (from numarray)
            #Note, if you are allowing directPhaseScreens, and wish not to bother creating an initial phase, use a dummy object in initPhs, with a copy method and a shape tuple.  This should then work saving you time and memory.
            self.phStructFn=None
            self.phStructFnIter=0
            self.coeffsZernike=None
            self.coeffsZernike2=None
            self.zernikeIters=0
            for key in self.parentList[0].keys():
                if type(initPhs[key]) in [type(""),type(0)]:#this is the idstr for the initial phasescreen...
                    #First, see if we can copy from parent - if not, create ourselves (this relies on the seed being a constant, not the current time).
                    try:
                        self.phaseScreens[key]=self.parentList[0][key].screen.copy()
                        print "infAtmos: Copied initial screen from parent %s"%str(key)
                    except:
                        print "infAtmos: cannot copy parent screen %s - generating directly."%str(key)
                        self.phaseScreens[key]=numpy.array(infScrn.computeInitialScreen(self.config,idstr=str(key)))
                else:
                    self.phaseScreens[key]=numpy.array(initPhs[key].copy())
                if self.phaseScreens[key].dtype.char!=self.dataType:
                    raise Exception("infAtmos and infScrn should use the same dataType value %s %s"%(self.phaseScreens[key].dtype.char,self.dataType))
            #Now get the layer offsets which are used to determine where in a layer to start looking for this source.  Basically, you can use sourceTheta and sourcePhi and sourceAlt to determine the distance in m from the on axis location.  Convert this into pixels, and you may get a negative number (ie if it is to the left of the onaxis position) depending on which source and which layer.  The maximum possible negative number is then equal to this layerOffset, which makes sure that all indexes into the array are positive.  This is best computed in the parameter file, as it requires knowledge of all sources, and of the size of each phase screen.   
            if self.atmosGeom==None:
                self.layerXOffset=self.config.getVal("layerXOffset")#dictionary depreciated
                self.layerYOffset=self.config.getVal("layerYOffset")#dictionary depreciated
            else:
                self.layerXOffset=self.atmosGeom.getLayerXOffset()
                self.layerYOffset=self.atmosGeom.getLayerYOffset()
            self.scrnXPxls={}
            self.scrnYPxls={}
            self.interpPhs=numpy.zeros((self.npup,self.npup),self.dataType)#resource sharing
            for key in self.parentList[0].keys():
                self.scrnXPxls[key]=self.phaseScreens[key].shape[1]#agbc from 0
                self.scrnYPxls[key]=self.phaseScreens[key].shape[0]#agbc from 1
            self.colAdd={}
            self.rowAdd={}
            self.maxColAdd={}
            self.maxRowAdd={}
            self.newCols={}
            self.newRows={}
            #create objects for working out the interpolation required for each new iteration...
            for key in self.parentList[0].keys():
                self.colAdd[key]=-self.vWind[key]*numpy.cos(self.windDirection[key]*numpy.pi/180)/self.pixScale*self.tstep#number of pixels to step each iteration (as float).  If <0, add at end.
                self.rowAdd[key]=-self.vWind[key]*numpy.sin(self.windDirection[key]*numpy.pi/180)/self.pixScale*self.tstep#number of pixels to step each iteration (as float).  If <0, add at end
                if self.rowAdd[key]==-0.0:
                    self.rowAdd[key]=0.
                if self.colAdd[key]==-0.0:
                    self.colAdd[key]=0.
                self.newCols[key]=util.getNewCols.getNewCols(numpy.fabs(self.colAdd[key]))
                self.newRows[key]=util.getNewCols.getNewCols(numpy.fabs(self.rowAdd[key]))
                self.maxColAdd[key]=int(numpy.ceil(numpy.fabs(self.colAdd[key])))
                self.maxRowAdd[key]=int(numpy.ceil(numpy.fabs(self.rowAdd[key])))
            # if colAdd<zero, we add new columns on the right of the array.
            # If rowAdd<zero, we add new rows at the bottom of the array (note, that if plotting in Gist, this is the top of the array).
            self.inputData={}#hold the new phase cols/rows.
            self.rowInput={}
            self.colInput={}
            self.interpPosCol={}
            self.interpPosRow={}
            self.nremCol={}
            self.nremRow={}
            self.naddCol={}
            self.naddRow={}
            for key in self.parentList[0].keys():
                self.inputData[key]=None
                self.rowInput[key]=None
                self.colInput[key]=None
            for i in xrange(len(self.idstr)):
                idstr=self.idstr[i]
                parent=self.parentList[i]
                self.initialise(parent,idstr)
    def newParent(self,parent,idstr=None):
        raise Exception("infAtmos - not yet able to accept new parent... (needs some extra coding)")
    def initialise(self,parentDict,idstr):
        """note, parent should be a dictionary here..."""
        if type(parentDict)!=type({}):#expecting a dictionary of parents
            parentDict={1:parentDict}
            print "WARNING - infAtmos: I was expecting a dictionary for parent here..."
        if self.parentKeyList==None:
            self.parentKeyList=parentDict.keys()
            self.nLayers=len(self.parentKeyList)
            self.parent=parentDict#store parent globally, since must all be same for every resource sharer.
        else:
            for key in parentDict.keys():
                if key not in self.parentKeyList:
                    raise Exception("infAtmos - all resource sharing objects must have same parents... (though whether they use them can be selected in the param file).")
        this=base.aobase.resourceSharer(parentDict,self.config,idstr,self.moduleName)
        self.thisObjList.append(this)
        
        if self.atmosGeom==None:
            sourceTheta=this.config.getVal("sourceTheta")#depreciated 
            sourceAlt=this.config.getVal("sourceAlt")#depreciated
            sourcePhi=this.config.getVal("sourcePhi")#depreciated
            sourceLam=this.config.getVal("sourceLam")
            zenith=0.
        else:
            sourceTheta=self.atmosGeom.sourceTheta(idstr)
            sourceAlt=self.atmosGeom.sourceAlt(idstr)
            sourcePhi=self.atmosGeom.sourcePhi(idstr)
            sourceLam=self.atmosGeom.sourceLambda(idstr)
            zenith=self.atmosGeom.zenith
            if sourceLam==None:
                sourceLam=self.config.getVal("sourceLam")
                print "Warning - atmosGeom sourceLam==None, using %g"%sourceLam
        layerList=this.config.getVal("atmosLayerList",default=None,raiseerror=0)
        intrinsicPhase=this.config.getVal("intrinsicPhase",default=None,raiseerror=0)
        storePupilLayers=this.config.getVal("storePupilLayers",default=0)
        computeUplinkTT=this.config.getVal("computeUplinkTT",default=0)
        launchApDiam=this.config.getVal("launchApDiam",default=0)
        this.atmosObj=util.atmos.atmos(parentDict,sourceAlt,sourceLam,sourceTheta,sourcePhi,self.npup,self.pupil,self.colAdd,self.rowAdd,self.layerAltitude,self.phaseScreens,self.scrnScale,self.layerXOffset,self.layerYOffset,layerList,zenith,intrinsicPhase=intrinsicPhase,storePupilLayers=storePupilLayers,computeUplinkTT=computeUplinkTT,launchApDiam=launchApDiam,ntel=self.ntel,telDiam=self.telDiam)

    def finalInitialisation(self):
        """since all memories are the same size (npup can't change...), its
        easy here..."""
        if self.doneFinalInit:
            return
        self.doneFinalInit=1
        for this in self.thisObjList:
            this.atmosObj.initMem(self.outputData,self.interpPhs)

    def generateNext(self,msg=None):
        """
        This function is called when it is okay to produce the next iteration
        of the simulation.
        Not expecting any msgs.
        """
        t1=time.time()
        if self.debug:
            print "infAtmos: GenerateNext (debug=%s)"%str(self.debug)
        if self.generate==1:
            if self.newDataWaiting:
                nvalid=0
                if self.currentIdObjCnt==0:
                    #if its the first object, we assume all screens are ready
                    #and construct them.   
                    for key in self.parent.keys():
                        if self.parent[key].dataValid==1:
                            if self.inputData[key] is not self.parent[key].outputData:
                                #create the arrays...
                                print "Allocating infAtmos arrays (hopefully this only happens during first iteration...)"
                                self.inputData[key]=self.parent[key].outputData
                                if self.inputData[key].dtype.char!=self.dataType:
                                    raise Exception("infAtmos and infScrn dataTypes must be same")
                                if self.parentSendWholeScreen==0:
                                    self.colInput[key]=self.inputData[key][:self.maxColAdd[key]*self.scrnYPxls[key]]#cmod.utils.arrayFromArray(self.inputData[key],(self.maxColAdd[key],self.scrnYPxls[key]),self.dataType)
                                    self.colInput[key].shape=(self.maxColAdd[key],self.scrnYPxls[key])
                                    self.rowInput[key]=self.inputData[key][self.maxColAdd[key]*self.scrnYPxls[key]:]#cmod.utils.arrayFromArray(self.inputData[key][self.maxColAdd[key]*self.scrnYPxls[key]:,],(self.maxRowAdd[key],self.scrnXPxls[key]),self.dataType)
                                    self.rowInput[key].shape=(self.maxRowAdd[key],self.scrnXPxls[key])
                            nvalid+=1
                    if nvalid==self.nLayers:
                        self.dataValid=1
                        self.makeLayers()#this is the first of resource sharers
                    elif nvalid>0:
                        print "ERROR: infAtmos - received wrong number of phase screens"
                        self.dataValid=0
                    else:
                        print "infAtmos: Waiting for data from infScrn, but not valid"
                        self.dataValid=0
                if self.dataValid:
                    if self.control["cal_source"]:#calibration source
                        self.outputData[:,]=0.
                    else:
                        self.thisObjList[self.currentIdObjCnt].atmosObj.createPupilPhs(self.phaseScreens,self.interpPosCol,self.interpPosRow,self.control)
                        if self.control["profilePhase"]:#compute phase covariance and profile.  This won't work if resource sharing.
                            #But doesn't matter, because its only really for testing anyway.
                            self.zernikeVariance(self.outputData,forDisplay=0)
                            self.phaseStructFunc(self.outputData,forDisplay=0)
                            
            else:#no new data ready
                self.dataValid=0
        else:
            self.dataValid=0
        if self.debug:
            print "infAtmos: done atmos generateNext (debug=%s)"%str(self.debug)
        self.generateNextTime=time.time()-t1


    def makeLayers(self):
        for key in self.parent.keys():#for each atmosphere layer...
            if self.directPhaseScreen and hasattr(self.parent[key],"screen"):
                #use the exact same copy from the infScrn object.
                #Note, this only works if it is in the same python process.
                self.phaseScreens[key]=self.parent[key].screen
            else:#construct the layer.
                self.makeLayer(key)
        
    def makeLayer(self,key):
        """make the phasescreen from the newly added cols/rows passed from infScrn object"""
        #first, add the new phase values into the phase array.
        nremCol,naddCol,self.interpPosCol[key]=self.newCols[key].next()
        nremRow,naddRow,self.interpPosRow[key]=self.newRows[key].next()
        self.nremCol[key]=nremCol
        self.naddCol[key]=naddCol
        self.nremRow[key]=nremRow
        self.naddRow[key]=naddRow
        if self.parentSendWholeScreen:
            self.phaseScreens[key]=self.parent[key].outputData
        else:
            shape=self.phaseScreens[key].shape
            #here, we add newly created phase from infScrn to the existing phasescreen stored here.
            if self.colAdd[key]<0:#add at end
                if self.rowAdd[key]<0:#add at end
                    self.phaseScreens[key][:shape[0]-naddRow,:shape[1]-naddCol]=self.phaseScreens[key][naddRow:,naddCol:,]
                    self.phaseScreens[key][:,shape[1]-naddCol:,]=numpy.transpose(self.colInput[key][self.maxColAdd[key]-naddCol:,])
                    self.phaseScreens[key][shape[0]-naddRow:,]=self.rowInput[key][self.maxRowAdd[key]-naddRow:,]
                else:#add at start
                    #interpPosRow=1-interpPosRow
                    #phsShiftX=1-phsShiftX
                    self.phaseScreens[key][naddRow:,:shape[1]-naddCol]=self.phaseScreens[key][:shape[0]-naddRow,naddCol:,].copy()
                    self.phaseScreens[key][:,shape[1]-naddCol:,]=numpy.transpose(self.colInput[key][self.maxColAdd[key]-naddCol:,])
                    self.phaseScreens[key][:naddRow]=self.rowInput[key][:naddRow]
            else:#add at start
                #interpPosCol=1-interpPosCol
                #phsShiftY=1-phsShiftY
                if self.rowAdd[key]<0:
                    self.phaseScreens[key][:shape[0]-naddRow,naddCol:,]=self.phaseScreens[key][naddRow:,:shape[1]-naddCol].copy()
                    self.phaseScreens[key][:,:naddCol]=numpy.transpose(self.colInput[key][self.maxColAdd[key]-naddCol:,])
                    self.phaseScreens[key][shape[0]-naddRow:,]=self.rowInput[key][self.maxRowAdd[key]-naddRow:,]
                else:
                    #interpPosRow=1-interpPosRow
                    #phsShiftX=1-phsShiftX
                    self.phaseScreens[key][naddRow:,naddCol:,]=self.phaseScreens[key][:shape[0]-naddRow,:shape[1]-naddCol].copy()
                    self.phaseScreens[key][:,:naddCol]=numpy.transpose(self.colInput[key][self.maxColAdd[key]-naddCol:,])
                    self.phaseScreens[key][:naddRow]=self.rowInput[key][:naddRow]


    def zernikeVariance(self,phase,forDisplay=0):
        """computes the zernike coefficients and stores them.
        phase here includes the pupil (though possibly doesn't have to).
        """
        rt=None
        if type(phase)!=type(None):
            if type(self.zernObj)==type(None):
                print "Creating Zernike object - may take a while (jmax=%d)"%self.jmax
                self.zernObj=util.zernikeMod.Zernike(self.pupil.fn,self.jmax,computeInv=1)
            coeffsZernike=self.zernObj.giveZernikeExpansion(phase);
            if type(self.coeffsZernike)==type(None):
                self.coeffsZernike=coeffsZernike[1:,].copy()
                self.coeffsZernike2=coeffsZernike[1:,]*coeffsZernike[1:,]
            else:
                self.coeffsZernike+=coeffsZernike[1:,]
                self.coeffsZernike2+=coeffsZernike[1:,]*coeffsZernike[1:,]
            self.zernikeIters+=1
        if forDisplay:
            #m1=numpy.average(coeffsZernike[ley][1:,]*coeffsZernike[key][1:,],axis=1)
            #m2=numpy.average(coeffsZernike[key][1:,],axis=1)
            varZ=self.coeffsZernike2/self.zernikeIters-(self.coeffsZernike*self.coeffsZernike)/self.zernikeIters**2
            #varZ=(m1-m2*m2)#get the data
            tabJ=numpy.arange(varZ.shape[0])+2#x axis
            # these can be plot using
            # fma();plmk(varZ,tabJ,marker=1,msize=0.5);plg(varTh,tabJ,marks=0);limits(2,jmax,1e-4,1e-1);logxy(1,1);
            rt=varZ,tabJ
            print "This doesn't seem to work - don't know why"
        return rt
    def phaseStructFunc(self,phase,forDisplay=0):
        """get the phase structure function.
        phase here doesn't include the pupil, ie is square phase, no mask.
        This should be called for lots of iterations to build up phStructFn.
        Then the result can be displayed.
        """
        rt=None
        if type(phase)!=type(None):
            dpix=phase.shape[0]
            ref=phase[0,0]#phase[dpix/2,dpix/2]
            if type(self.phStructFn)==type(None):
                self.phStructFn=(phase-ref)**2.#
                self.phaseCovariance2=phase*ref
                self.phaseSum=phase.copy()
                self.phaseSumCnt=1
                self.phStructFnIter=0
            else:
                self.phStructFn+=(phase-ref)**2.
                self.phaseCovariance2+=phase*ref
                self.phaseSum+=phase
                self.phaseSumCnt+=1
                self.phaseCovariance=self.phaseCovariance2/self.phaseSumCnt-(self.phaseSum*self.phaseSum[0,0])/(self.phaseSumCnt*self.phaseSumCnt)
            self.phStructFnIter+=1
        if forDisplay:
            pr=numpy.array(util.dist.profilRadial(self.phStructFn/self.phStructFnIter));
            rr=pr[0,]*self.pixScale
            pr=pr[1,]
            rt=pr,rr
            print "This doesn't seem to work - don't know why - maybe because the pupil mask is inplace?"
            #window(1);logxy(0,0);fma();plg(pr,rr);l=limits()
        return rt
    
    def theoreticalZernikeVariance(self,L0=None,r0=None):
        """get theoretical zernike variance.  
        This is typically called by the GUI"""
        if L0==None:
            L0=self.L0
        if r0==None:
            r0=self.r0
        varTh=numpy.array(util.zernikeMod.matCovZernikeL0(self.jmax, self.telDiam/L0, self.telDiam/r0,diagOnly=1)[1:,])#phase covariance
        #varTh=numpy.diagonal(matCovZ);#diagonal is the phase variance.
        return varTh
    def theoreticalPhaseStructFunc(self,L0=None,r0=None,npup=None):
        """Computation of the theoretical phase structure function.
        Typically called by the GUI."""
        if L0==None:
            L0=self.L0
        if r0==None:
            r0=self.r0
        if npup==None:
            npup=self.ntel
        f0=1./L0;
        distMap=numpy.array(util.dist.profilRadial(None,radiusOnly=1,nfft=npup,dtype="d"))*self.pixScale
        coeff=0.17166*(r0*f0)**(-5./3)
        K=gamma(5./6)/(2**(1./6))#gamma function from scipy.special
        DPhiTh=(2*numpy.pi*distMap*f0)**(5./6)
        DPhiTh*=kv(5./6,2*numpy.pi*distMap*f0)#bessel function from scipy.special
        DPhiTh=coeff*(K-DPhiTh);DPhiTh[0]=0
        #plg(DPhiTh,distMap,color="red")
        return DPhiTh

    def createSingleLayerPhase(self,key):
        """Primarily for GUI... - gets the phase for this pupil at one height."""
        data=self.thisObjList[self.currentIdObjCnt].atmosObj.createSingleLayerPhs(self.phaseScreens,self.interpPosCol,self.interpPosRow,key,self.control)
        return data


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
        txt="""<plot title="Pupil phase screen%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
        if len(self.thisObjList)==1:
            txt+="""<plot title="profile phase%s?" cmd="feedback=1-%s.control['profilePhase'];%s.control['profilePhase']=feedback" ret="feedback" when="cmd" texttype="1" wintype="mainwindow">\nbutton=feedback\n</plot>"""%(id,objname,objname)
            txt+="""<plot title="phase structure function%s" cmd="pr,rr=%s.phaseStructFunc(None,1);data=numpy.array((rr,pr,%s.theoreticalPhaseStructFunc()))" ret="data" type="pylab" when="rpt" dim="1"/>"""%(id,objname,objname)
            txt+="""<plot title="zernike variance %s" cmd="varZ,tabJ=%s.zernikeVariance(None,1);data=numpy.array((tabJ.astype('f'),varZ,%s.theoreticalZernikeVariance()))" ret="data" type="pylab" when="rpt" dim="1"/>"""%(id,objname,objname)
            for key in self.parent.keys():
                txt+="""<plot title="layer phase %s%s" cmd="data=%s.createSingleLayerPhase('%s')" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(key,id,objname,key)

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
        #paramList.append(base.dataType.dataType(description="degRad",typ="eval",val="2*numpy.pi/360.",comment="degrees to radians."))        
        #paramList.append(base.dataType.dataType(description="arcsecRad",typ="eval",val="2*numpy.pi/360/3600",comment="arcsec to radians."))        
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="ntel",typ="eval",val="this.globals.npup",comment="Pixels for telescope"))
        #paramList.append(base.dataType.dataType(description="telSec",typ="f",val="8.",comment="TODO: Telescope secondary diameter (m)"))
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(this.globals.npup,this.globals.ntel/2,this.globals.ntel/2*this.globals.telSec/this.globals.telDiam,this.globals.wfs_nsubx,this.globals.wfs_minarea)",comment="Telescope pupil"))
        #paramList.append(base.dataType.dataType(description="windDirection",typ="eval",val="{'0m':0.,'2000m':60.}",comment="TODO: Wind direction (degrees, going from -180 to 180) - dict with keys equal to layer heights."))
        #paramList.append(base.dataType.dataType(description="vWind",typ="eval",val="{'0m':10.,'2000m':13.}",comment="TODO: Wind velocity (m/s) - dict with keys equal to layer heights."))
        #paramList.append(base.dataType.dataType(description="altitude",typ="eval",val="{'0m':0.,'2000m':2000.}",comment="TODO: Layer altitudes - dict with keys equal to layer heights."))
        #paramList.append(base.dataType.dataType(description="sourceThetaDict",typ="eval",val="{'onaxis':0.,'offaxis':600.}",comment="TODO: Source positions (theta)"))
        #paramList.append(base.dataType.dataType(description="sourcePhiDict",typ="eval",val="{'onaxis':0.,'offaxis':30.}",comment="TODO: Source positions (phi)"))
        #paramList.append(base.dataType.dataType(description="sourceTheta",typ="eval",val="this.infAtmos.sourceThetaDict['onaxis']",comment="TODO: Source positions (theta)"))
        #paramList.append(base.dataType.dataType(description="sourcePhi",typ="eval",val="this.infAtmos.sourcePhiDict['onaxis']",comment="TODO: Source positions (phi)"))
        #paramList.append(base.dataType.dataType(description="sourceLamDict",typ="eval",val="{'onaxis':1650.,'offaxis':1650.}",comment="TODO: Source wavelength"))
        #paramList.append(base.dataType.dataType(description="sourceLam",typ="eval",val="this.infAtmos/sourceLamDict['onaxis']",comment="TODO: Source wavelength"))
        #paramList.append(base.dataType.dataType(description="sourceAltDict",typ="eval",val="{'onaxis':-1.,'offaxis':-1.}",comment="TODO: Source height (positive for LGS)"))
        #paramList.append(base.dataType.dataType(description="sourceAlt",typ="eval",val="this.infAtmos.sourceAltDict['onaxis']",comment="TODO: Source height (positive for LGS)"))
        paramList.append(base.dataType.dataType(description="tstep",typ="f",val="0.005",comment="TODO: timestep."))        
        #paramList.append(base.dataType.dataType(description="layerOffset",typ="code",val="from science.infAtmos import calcLayerOffset;layerOffset=calcLayerOffset(this.infScrn.scrnSize,this.infAtmos.sourceThetaDict,this.infAtmos.sourcePhiDict,this.infAtmos.altitude,this.globals.npup,this.globals.ntel,this.globals.telDiam)",comment="Layer offsets"))
        #paramList.append(base.dataType.dataType(description="layerXOffset",typ="eval",val="this.infAtmos.layerOffset['onaxis']",comment="TODO: Layer offset"))
        #paramList.append(base.dataType.dataType(description="layerYOffset",typ="eval",val="this.infAtmos.layerOffset['onaxis']",comment="TODO: Layer offset"))
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))
        return paramList
