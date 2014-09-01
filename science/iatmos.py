#"""Module iatmos - takes one or more iscrn outputs and combines together,
#does interpolation etc and creates the telescope pupil phase.
#Using the newer, memory efficient iscrn module.
#"""

import base.aobase,util.getNewCols
import numpy
import time,types
import iscrn
import util.atmos
import util.dist,util.zernikeMod
from scipy.special import gamma,kv
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


class iatmos(base.aobase.aobase):
    """Create an iatmos object.  This object can take several iscrn
    objects as parents, and returns the pupil phase for a given source

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
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Initialise the object.  The parent object here should be a dictionary with values for each atmospheric layer, ie instances of iscrn.  """
        if type(parent)!=type({}):
            parent={"L0":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.parentDict=parent
        self.scrnDataType="d"#self.config.getVal("dataType")
        self.outDataType="f"
        self.npup=self.config.getVal("npup")#pupil size (pixels)
        if forGUISetup==1:
            self.outputData=[(self.npup,self.npup),self.outDataType]
        else:#setup for real
            self.timing=self.config.getVal("timing",default=0)
            self.interpolationNthreads=self.config.getVal("interpolationNthreads",default=0)#a tuple of (nthreads,nblockx,nblocky)
            self.parentSendWholeScreen=0
            self.doneFinalInit=0
            self.sentPlotsCnt=0
            self.layerListDict={}
            self.atmosGeom=self.config.getVal("atmosGeom")
            self.directPhaseScreen=self.config.getVal("directPhaseScreen",default=1)#are we allowed to access parent.screen if parent is an iscrn object?  
            self.telDiam=self.config.getVal("telDiam")
            self.ntel=self.config.getVal("ntel")#tel diameter in pxls.  
            self.pupil=self.config.getVal("pupil")
            #self.telSec=self.config.getVal("telSec")
            self.scrnScale=self.telDiam/float(self.ntel)#self.config.getVal("scrnScale")
            scrnScale=self.scrnScale
            self.pixScale=self.telDiam/float(self.npup)
            self.jmax=self.config.getVal("jmax",default=55)#jmax here is by default 55.  A better value might be npup*0.75 ^2/2, about the most you would hope.  If npup is large, this will take a long time, so you should probably use less for jmax.  Only used for zernike variance (from GUI).
            self.zernObj=None
            #The number of parent objects determines the number of layers.
            #self.nLayers=len(self.parentList[0].keys())
            self.windDirection=self.atmosGeom.windDirections()
            self.vWind=self.atmosGeom.windSpeeds()
            self.layerAltitude=self.atmosGeom.layerAltitudes()#these are scaled by zenith.
            self.r0=self.atmosGeom.r0
            self.L0=self.atmosGeom.l0
    
            self.control={"cal_source":0,"profilePhase":0,"fullPupil":0,"removePiston":1}#full pupil is used as a flag for phase profiling (xinterp_recon etc) if want to return the whole pupil.
            self.tstep=self.config.getVal("tstep")
            #self.niters=0
            self.outputData=numpy.zeros((self.npup,self.npup),self.outDataType)#resource sharing
            #self.scrnXPxls=self.config.getVal("scrnXPxls")#will depend on which screen it is, so this is no good!
            #self.scrnYPxls=self.config.getVal("scrnYPxls")
            #Get the initial phase screen.
            self.phaseScreens={}
            self.nLayers=len(parent.keys())

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
            for pkey in self.parent.keys():
                #each parent can have 1 or many layers.  So, get this into a dictionary.
                self.layerListDict[pkey]=self.config.getVal("layerList",{},searchOrder=["iscrn_%s","iscrn","globals"]).get(pkey,[pkey])


                for key in self.layerListDict[pkey]:
                    #First, see if we can copy from parent - if not, create ourselves (this relies on the seed being a constant, not the current time).
                    try:
                        self.phaseScreens[key]=self.parent[pkey].thisObjDict[key].screen#[xxx].copy()
                        print "iatmos: Copied initial screen from parent %s[%d]"%(str(pkey),str(key))
                        print "TODO: iatmos - is this copy of initial screen needed?"
                    except:
                        print "iatmos: cannot copy parent screen %s - generating directly."%str(key)
                        self.phaseScreens[key]=numpy.array(iscrn.computeInitialScreen(self.config,idstr=str(key)))
                    if self.phaseScreens[key].dtype.char!=self.scrnDataType:
                        raise Exception("iatmos and iscrn should use the same dataType value %s %s"%(self.phaseScreens[key].dtype.char,self.scrnDataType))
                    ps=self.phaseScreens[key]
                    if ps.shape!=(atmosGeom.getScrnYPxls(key,rotateDirections=1),atmosGeom.getScrnXPxls(key,rotateDirections=1)):
                        raise Exception("Phase screen size unexpected in iatmos: %s -> %s"%(str(ps.shape),str((atmosGeom.getScrnYPxls(key,rotateDirections=1),atmosGeom.getScrnXPxls(key,rotateDirections=1)))))
                    if self.ygradient!=None:
                        self.ygradient[key]=numpy.empty(ps.shape,self.scrnDataType)
                        self.ygradient[key][1:-1]=(ps[2:]-ps[:-2])*0.5
                        self.ygradient[key][0]=ps[1]-ps[0]
                        self.ygradient[key][-1]=ps[-1]-ps[-2]

            #Now get the layer offsets which are used to determine where in a layer to start looking for this source.  Basically, you can use sourceTheta and sourcePhi and sourceAlt to determine the distance in m from the on axis location.  Convert this into pixels, and you may get a negative number (ie if it is to the left of the onaxis position) depending on which source and which layer.  The maximum possible negative number is then equal to this layerOffset, which makes sure that all indexes into the array are positive.  This is best computed in the parameter file, as it requires knowledge of all sources, and of the size of each phase screen.   
            self.layerXOffset=self.atmosGeom.getLayerXOffset(rotateDirections=1)
            self.layerYOffset=self.atmosGeom.getLayerYOffset(rotateDirections=1)
            self.scrnXPxls={}
            self.scrnYPxls={}
            #self.interpPhs=numpy.zeros((self.npup,self.npup),self.dataType)#resource sharing
            for pkey in self.parent.keys():
                for key in self.layerListDict[pkey]:
                    self.scrnXPxls[key]=self.phaseScreens[key].shape[1]#agbc from 0
                    self.scrnYPxls[key]=self.phaseScreens[key].shape[0]#agbc from 1
            self.rowAdd={}
            self.maxRowAdd={}
            self.newRows={}
            self.insertPos={}
            #create objects for working out the interpolation required for each new iteration...
            for pkey in self.parent.keys():
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
            for key in self.parent.keys():
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
            print "WARNING - iatmos: I was expecting a dictionary for parent here..."
        if not set(parentDict)==set(self.parentDict):#check they have teh same parents (main and resource sharer).
            raise Exception("iatmos - all resource sharing objects must have same parents... (though whether they use them can be selected in the param file")
        this=base.aobase.resourceSharer(parentDict,self.config,idstr,self.moduleName)
        self.thisObjList.append(this)
        
        sourceTheta=self.atmosGeom.sourceTheta(idstr)
        sourceAlt=self.atmosGeom.sourceAlt(idstr)
        sourcePhi=self.atmosGeom.sourcePhi(idstr)
        sourceLam=self.atmosGeom.sourceLambda(idstr)
        zenith=self.atmosGeom.zenith
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
        this.atmosObj=util.atmos.iatmos(sourceAlt,sourceLam,sourceTheta,sourcePhi,self.npup,self.pupil,self.rowAdd,self.layerAltitude,self.windDirection,self.phaseScreens,self.ygradient,self.scrnScale,self.layerXOffset,self.layerYOffset,layerList,zenith,intrinsicPhase=intrinsicPhase,storePupilLayers=storePupilLayers,computeUplinkTT=computeUplinkTT,launchApDiam=launchApDiam,ntel=self.ntel,telDiam=self.telDiam,interpolationNthreads=self.interpolationNthreads,outputData=self.outputData)

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
        if self.debug:
            print "iatmos: GenerateNext (debug=%s)"%str(self.debug)
        if self.generate==1:
            if self.newDataWaiting:
                nvalid=0
                if self.currentIdObjCnt==0:
                    #if its the first object, we assume all screens are ready
                    #and construct them.   
                    for pkey in self.parent.keys():
                        if self.parent[pkey].dataValid==1:
                            if self.inputData[pkey] is not self.parent[pkey].outputData:
                                #create the arrays...
                                #print "Allocating iatmos arrays (hopefully this only happens during first iteration...)"
                                self.inputData[pkey]=self.parent[pkey].outputData
                                if self.inputData[pkey].dtype.char!=self.scrnDataType:
                                    raise Exception("iatmos and iscrn dataTypes must be same")
                                if self.parentSendWholeScreen==0:
                                    pos=0
                                    for key in self.layerListDict[pkey]:
                                        self.rowInput[key]=self.inputData[pkey][pos:pos+self.maxRowAdd[key]*self.scrnXPxls[key]]
                                        pos+=self.maxRowAdd[key]*self.scrnXPxls[key]
                                        self.rowInput[key].shape=(self.maxRowAdd[key],self.scrnXPxls[key])
                            nvalid+=1
                    if nvalid==len(self.parent.keys()):#self.nLayers:
                        self.dataValid=1
                        self.makeLayers()#this is the first of resource sharers
                    elif nvalid>0:
                        print "ERROR: iatmos - received wrong number (%d/%d) of phase screens"%(self.nvalid,len(self.parent.keys()))#self.nLayers)
                        self.dataValid=0
                    else:
                        print "iatmos: Waiting for data from iscrn, but not valid"
                        self.dataValid=0
                if self.dataValid:
                    if self.control["cal_source"]:#calibration source
                        self.outputData[:,]=0.
                    else:
                        self.thisObjList[self.currentIdObjCnt].atmosObj.createPupilPhs(self.interpPosRow,self.insertPos,self.control)
                        if self.control["profilePhase"]:#compute phase covariance and profile.  This won't work if resource sharing.
                            #But doesn't matter, because its only really for testing anyway.
                            self.zernikeVariance(self.outputData,forDisplay=0)
                            self.phaseStructFunc(self.outputData,forDisplay=0)
            else:#no new data ready
                self.dataValid=0
        else:
            self.dataValid=0
        if self.debug:
            print "iatmos: done atmos generateNext (debug=%s)"%str(self.debug)
        self.generateNextTime=time.time()-t1


    def makeLayers(self):
        for pkey in self.parent.keys():#for each atmosphere layer...
            if self.directPhaseScreen and hasattr(self.parent[pkey],"thisObjDict"):
                #we can just use the screens from the parents (though we'll still have to create the gradients).
                for key in self.layerListDict[pkey]:
                    
                    #use the exact same copy from the iscrn object.
                    #Note, this only works if it is in the same python process.
                    self.phaseScreens[key]=self.parent[pkey].thisObjDict[key].screen
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
            self.phaseScreens[key]=self.parent[key].outputData
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
        txt="""<plot title="Pupil phase screen%s (last resource sharer)" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
        if len(self.thisObjList)==1:
            txt+="""<plot title="profile phase%s?" cmd="feedback=1-%s.control['profilePhase'];%s.control['profilePhase']=feedback" ret="feedback" when="cmd" texttype="1" wintype="mainwindow">\nbutton=feedback\n</plot>"""%(id,objname,objname)
            txt+="""<plot title="phase structure function%s" cmd="pr,rr=%s.phaseStructFunc(None,1);data=numpy.array((rr,pr,%s.theoreticalPhaseStructFunc()))" ret="data" type="pylab" when="rpt" dim="1"/>"""%(id,objname,objname)
            txt+="""<plot title="zernike variance %s" cmd="varZ,tabJ=%s.zernikeVariance(None,1);data=numpy.array((tabJ.astype('f'),varZ,%s.theoreticalZernikeVariance()))" ret="data" type="pylab" when="rpt" dim="1"/>"""%(id,objname,objname)
            for pkey in self.parent.keys():
                for key in self.layerListDict[pkey]:
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
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="ntel",typ="eval",val="this.globals.npup",comment="Pixels for telescope"))
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(this.globals.npup,this.globals.ntel/2,this.globals.ntel/2*this.globals.telSec/this.globals.telDiam,this.globals.wfs_nsubx,this.globals.wfs_minarea)",comment="Telescope pupil"))
        paramList.append(base.dataType.dataType(description="tstep",typ="f",val="0.005",comment="TODO: timestep."))        
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by iscrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (iatmos), theta, phi, alt, nsubx or None)."))
        return paramList
