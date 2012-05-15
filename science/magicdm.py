import numpy
import util.zernikeMod
import util.FITS
import base.aobase
class dm(base.aobase.aobase):
    """
    Zernike DM fitter:
    Takes the input phase, decomposes it to a zernike basis, and then reconstructs using this basis.
    This DM does not need a reconstructor.
    Add DM phases to input phase array and output.
    Note, this probably won't work very well when resource sharing since the dm is shaped from the first atmos input only.
    If needed, this could easily be updated to include other modes, eg actuator pokes/annular zernikes etc.
    Also can compute zernike covariances.
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        if forGUISetup==1:
            npup=self.config.getVal("npup")
            self.outputData=[(npup,npup),numpy.float64]
        else: # set up for simulation.
            self.npup=self.config.getVal("npup")
            self.atmosGeom=self.config.getVal("atmosGeom",default=None,raiseerror=0)
            self.dmObj=self.config.getVal("dmObj",default=None,raiseerror=0)
	    self.dmpup=self.dmObj.calcdmpup(self.idstr[0])#number of pixels to store the phase. May be >npup if not ground conjugate.
	    self.conjHeight=self.dmObj.getHeight(self.idstr[0])
            self.pupil=self.config.getVal("pupil")

            self.nmodes=config.getVal("nmodes")#number of zernikes to use.
            self.covModes=config.getVal("covModes",default=self.nmodes)#number of modes to use when computing covariances.
            self.coeff=numpy.zeros((self.nmodes),numpy.float64)
            self.covCoeff=numpy.zeros((self.covModes),numpy.float64)
            self.monteCovMat=None
            self.montePhaseCovMat=None#this one should have same scaling as Noll in xinterp_recon
            self.dmphs=numpy.zeros((self.dmpup,self.dmpup),numpy.float32)            # DM figure      
            #self.gamma=config.getVal("gamma")#gain for zernikes.  Can be a float or an array length nmodes.
            #There are 2 ways of getting zernikes - the original, using RWW way, which returns circular zernikes, or FA way, which allows a pupil mask to be specified (and so zernikes aren't always orthogonal)
            self.sourceLam=self.dmObj.getDM(self.idstr[0]).reconLam#atmosGeom.sourceLambda(self.idstr[0])#self.config.getVal("sourceLam")#wavelength for which the mirror will be shaped.
            #self.wfslam=self.config.getVal("wfslam")#reconstructor wavelength
	    self.zoffset=self.config.getVal("zoffset",default=None,raiseerror=0)#zernike offsets that can be applied.  An array of length nmodes, ie numpy.zeros((self.nmodes,),numpy.float64)

	    self.reconData=None
            #self.wavelengthRatio=self.wfslam/self.sourcelam
            self.telDiam=self.config.getVal("telDiam")
            self.outputData=numpy.zeros((self.npup,self.npup),numpy.float64)

            # Make Zernike fns over DM, and
            if max(self.nmodes,self.covModes)>0:
                self.zern=util.zernikeMod.Zernike(self.pupil,max(self.nmodes,self.covModes),computeInv=0).zern
                util.zernikeMod.normalise(self.zern)#normalise to orthonormal (ie numpy.sum(zern[i]*zern[i])==1).
            self.control={"dm_update":1,"computeCovariance":0}

            for i in xrange(len(self.idstr)):
                self.initialise(self.parentList[i],self.idstr[i])


    def initialise(self,parentDict,idstr):
        """note parent may be a dictionary here (atmos and recon)"""
        if type(parentDict)!=type({}):
            parentDict={"atmos":parentDict}
        this=base.aobase.resourceSharer(parentDict,self.config,idstr,self.moduleName)
        self.thisObjList.append(this)
        #for a source with theta and a conjugate height c, the
        #separation at this height will be c tan(theta).  As a
        #fraction of the phasescreen width, this will be c
        #tan(theta)/telDiam.  So, in pixels, this will be c
        #tan(theta)/telDiam * npup.  Multiply this by
        #cos/sin(phi), and you then have the correct offsets.
        this.xoff=(self.dmpup-self.npup)/2#central in the case of
        this.yoff=(self.dmpup-self.npup)/2#ground conjugated.
        this.sourceID=self.dmObj.getSourceID(idstr)
        this.sourceAlt=self.atmosGeom.sourceAlt(this.sourceID)

        if self.conjHeight!=0:
            print "todo: check that signs are right - and that xoff and yoff are the right way round... (zdm - conjugated height)"
            this.sourceTheta=self.atmosGeom.sourceTheta(this.sourceID)*numpy.pi/180/3600.
            this.sourcePhi=self.atmosGeom.sourcePhi(this.sourceID)*numpy.pi/180
            r=self.conjHeight*numpy.tan(this.sourceTheta)/self.telDiam*self.npup
            this.xoff+=int(r*numpy.cos(this.sourcePhi)+0.5)#round
            this.yoff+=int(r*numpy.sin(this.sourcePhi)+0.5)#correctly
            #but don't bother with any interpolation here...
        # select the correct part of the dm phs.
        if this.xoff<0 or this.yoff<0 or this.xoff+self.npup>self.dmpup or this.yoff+self.npup>self.dmpup:
            raise Exception("DM pupil not large enough to hold all the conjugated sources... %d %d %d %d"%(this.xoff,this.yoff,self.npup,self.dmpup))
        this.selectedDmPhs=self.dmphs[this.yoff:this.yoff+self.npup,this.xoff:this.xoff+self.npup]
        this.sourceLam=self.atmosGeom.sourceLambda(this.sourceID)#this.config.getVal("sourceLam")          # Wavelength for this optical path
        #this.wfsLam=this.config.getVal("wfslam",default=this.sourceLam)
        #this.wavelengthRatio=this.wfsLam/this.sourceLam#wfs_lam/lam
        this.wavelengthAdjustor=self.sourceLam/this.sourceLam#this.wavelengthRatio/self.wavelengthRatio#the mirror will be shaped as for self.wavelengthRatio... if this.sourceLam is longer than self.sourceLam, radians of phase P-V will be smaller, so less change needed in wavelengthAdjustor.


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
        self.thisObjList[indx].parent=parent
        
    def generateNext(self,ms=None):
        """DM main loop:  Reflect from DM (add DM figure to phase)"""
        this=self.thisObjList[self.currentIdObjCnt]
        if not this.parent.has_key("atmos"):
            print "magicdm object assigning parents automatically"
            keylist=this.parent.keys()
            for key in keylist:
                if this.parent[key].outputData.shape==(self.npup,self.npup):#the right shape for an atmos object
                    print "magicdm parent object %s becoming atmos with output shape %s"%(str(key),str(this.parent[key].outputData.shape))
                    this.parent["atmos"]=this.parent[key]
                else:
                    print "magicdm parent object %s being ignored withh output shape %s"%(str(key),str(this.parent[key].outputData.shape))
            for key in keylist:
                del(this.parent[key])
        if not this.parent.has_key("atmos"):#maybe just running for poking?
            print "magicdm object no parent atmos object found.  Assuming unperturbed phase in - will do nothing!"
            class dummyAtmosClass:
                dataValid=1
                outputData=0
            this.parent["atmos"]=dummyAtmosClass()
        if self.generate==1:
            if self.newDataWaiting:
                if this.parent["atmos"].dataValid==1:
                    self.outputData[:,]=this.parent["atmos"].outputData
                    self.dataValid=1
                else:
                    print "DM: waiting for phase data, but not valid"
                    self.dataValid=0
                if self.control["dm_update"]==1 and self.currentIdObjCnt==0:
                    self.update()
                    self.dataValid=1#update the output.
            if self.dataValid:
                self.selectedDmPhs=this.selectedDmPhs
                if this.wavelengthAdjustor==1:#dm is shaped for this wavelength...
                    self.outputData+=self.selectedDmPhs
                else:
                    self.outputData+=self.selectedDmPhs*this.wavelengthAdjustor
                self.outputData*=self.pupil.fn
                #now remove piston
                #phasesum=numpy.sum(self.outputData.ravel())
                #self.outputData-=phasesum/self.pupil.sum
                #self.outputData*=self.pupil.fn
                #self.outputData += (self.dmphs*(self.wavelengthRatio))# Calculate reflected phase
                if self.control["computeCovariance"]:
                    self.computeCovariance(self.outputData)
        else:
            self.dataValid=0
            
    def update(self):
        """DM figure update: Add zernike updates to
        mirror (if the loop is closed)"""
	self.dmphs[:]=0#numpy.zeros((self.npup,self.npup),numpy.float64)
        #zoffset is only used here (I think).  Should be set by GUI.
        #zpoke also used by reconstructor (but doesn't do anything there)
        #powers=self.gamma*self.reconData
        self.coeff[:]=0
        phs=self.thisObjList[0].parent["atmos"].outputData
	for j in range(self.nmodes):
            self.coeff[j]=-numpy.sum(phs*self.zern[j])#compute the coefficient for this zernike...
	    self.dmphs+=self.coeff[j]*self.zern[j]     # and add it to the mirror shape
	    
        if type(self.zoffset)!=type(None):
            for j in range(min(self.nmodes,self.zoffset.shape[0])):
                self.dmphs+=self.zoffset[j]*self.zern[j]#always add constant offsets or pokes
                
    def computeCovariance(self,phs):
        for i in xrange(self.covModes):
            self.covCoeff[i]=numpy.sum(phs*self.zern[i])
        if type(self.monteCovMat)==type(None):
            self.monteCovMat=numpy.zeros((self.covModes,self.covModes),numpy.float64)
            self.montePhaseCovMat=numpy.zeros((self.covModes,self.covModes),numpy.float64)
            self.navCovMat=0
        self.monteCovMat+=self.covCoeff[:,numpy.newaxis]*self.covCoeff[numpy.newaxis,:]
        self.navCovMat+=1
        self.montePhaseCovMat[:]=self.monteCovMat/(self.navCovMat*self.pupil.sum)



    
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
            txt+="""<plot title="magicdm output%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="magicdm mirror%s" cmd="data=-%s.dmphs" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="magicdm coeffs%s" cmd="data=-%s.coeff" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="magicdm compute covariance%s" when="cmd" ret="feedback" texttype="1" wintype="mainwindow">\n<cmd>\nfeedback=1-%s.control['computeCovariance']\n%s.control['computeCovariance']=feedback\n</cmd>\nbutton=feedback\nif feedback==0:\n msg='Phase covariance finished computing'\nelse:\n msg='Computing phase covariance...'\nfeedback=msg</plot>\n"""%(id,objname,objname)
            txt+="""<plot title="phase covariance%s" cmd="data=%s.montePhaseCovMat" ret="data" type="pylab" when="rpt"/>\n"""%(id,objname)
        txt+="""<plot title="magicdm selected mirror%s" cmd="data=-%s.thisObjList[%d].selectedDmPhs" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt)
        self.sentPlotsCnt=(self.sentPlotsCnt+1)%len(self.thisObjList)
        return txt


