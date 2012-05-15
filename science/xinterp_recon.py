#import Numeric
import os
import numpy
import numpy.linalg
#from LinearAlgebra import generalized_inverse,singular_value_decomposition,inverse
#from plwf import *
import Scientific.MPI
import base.aobase
import util.FITS
import util.fitModes
import util.centroid
import cmod.phaseCov
import util.phaseCovariance
import util.fdpcg
import util.createPokeMx
import util.tel
import util.zernikeMod
try:
    import util.regularisation
except:
    print "WARNING - regularisation not supported - please update your scipy library"
from util.centroid import pxlToRadPistonTheoretical
from util.dm import MirrorSurface
#import numpy
import string

#print "Search for TODO in xinterp_recon.py"
class recon(base.aobase.aobase):
    """A reconstructor capable of several different types of reconstruction,
    which sends actuator values for a model xinetics DM to the children.
    
    The reconstruction types include a fourier domain PCG reconstructor, which
    is able to give about the same level of correction as a SOR reconstructor.
    It can also reconstruct MAP/zernikes or MVM.
    
    Currently classical AO only.  Use tomoRecon for others...
    """
    
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        global util
        if type(parent)!=type({}):
            parent={"cent":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        #if forGUISetup==1:
        #    nsubx=self.config.getVal("wfs_nsubx")
        #    print "xinterp_recon: Dimensions not yet known"
        #else: #set up the module for use...
        self.dataValid=1#data is valid even before the first iteration because we assume the mirror starts zerod
        self.pupil=self.config.getVal("pupil")

        self.atmosGeom=self.config.getVal("atmosGeom",default=None,raiseerror=0)
        self.dmObj=self.config.getVal("dmObj",default=None,raiseerror=0)
        if self.dmObj==None or self.atmosGeom==None:#depreciated mode
            print "xinterp_recon - dmObj or atmosGeom object not found in config - using depreciated mode."
            nact=self.nact=config.getVal("nAct")
            npup=self.npup=config.getVal("npup")
            nsubx=self.wfs_nsubx=config.getVal("wfs_nsubx")
            n=self.wfs_n=config.getVal("wfs_n")
            self.telDiam=config.getVal("telDiam")
            self.dmflag=self.config.getVal("dmflag")
            self.wfslambda=config.getVal("wfslam")
            self.interpType=self.config.getVal("dmInterpType")#type of interpolation done.
            self.actCoupling=self.config.getVal("actCoupling",raiseerror=0)
            self.actFlattening=self.config.getVal("actFlattening",raiseerror=0)
            self.r0=config.getVal("r0")#I think r0 is defined at 500nm.
            self.l0=config.getVal("l0")

        else:#use dmObj
            #Since xinterp_recon only works with classical AO, we here assume that there is only 1 DM, and get its idstr.
            self.r0=self.atmosGeom.r0
            self.l0=self.atmosGeom.l0
            self.dmidstr=self.config.getVal("dmidstr",default=self.idstr[0],raiseerror=0)
            if self.dmidstr==None:
                self.dmidstr=self.dmObj.dmInfoList[0].idlist[0][0]
                print "WARNING - using %s for dmidstr (xinterp_recon %s)"%(self.dmidstr,self.idstr[0])
            self.sourceid=self.config.getVal("sourceid",default=self.dmObj.getDM(self.dmidstr).idlist[0][1])
            #self.dmpup=self.dmObj.calcdmpup(dmidstr)
            self.dm=self.dmObj.getDM(self.dmidstr)
            if self.dm.zonalDM:
                nact=self.nact=self.dmObj.getnAct(self.dmidstr)
                self.actoffset=self.dmObj.getDM(self.dmidstr).actoffset
            else:
                nact=self.nact=self.dm.nact
                self.actoffset=None
            self.actCoupling=self.dm.actCoupling
            self.actFlattening=self.dm.actFlattening
            npup=self.npup=self.atmosGeom.npup
            nsubx=self.wfs_nsubx=self.atmosGeom.sourcensubx(self.sourceid)
            n=self.wfs_n=npup/nsubx
            self.telDiam=self.atmosGeom.telDiam
            self.dmflag=self.dmObj.computeDMPupil(self.dmidstr,centObscuration=self.pupil.r2,retPupil=0)[0]
            self.wfslambda=self.atmosGeom.sourceLambda(self.sourceid)
            if self.wfslambda==None:
                self.wfslambda=config.getVal("wfslam")
                print "source lambda not found in atmosGeom, using %g instead"%self.wfslambda
            self.interpType=self.dm.interpType
        self.nfft=config.getVal("wfs_nfft")
        self.nimg=self.config.getVal("wfs_nimg")
        self.rcond=config.getVal("rcond")
        self.gamma=config.getVal("gamma")#single value or array...
        self.removeAverageCent=self.config.getVal("removeAverageCent",default=0)#for GLAO...
        self.reconType=self.config.getVal("recontype")
        self.pokeval=config.getVal("pokeval",default=5.)#value to place on actuators when poking.  Needs to be eg 150 for zernDM dmModeType when used with zdm (large since the zernikes are normalised).
        self.monteNoll=None
        self.navNoll=0
        self.tiptilt=numpy.zeros(2,numpy.float64)# Closed loop tilt increments
        self.poke=0
        self.control={"poke":0,"close_dm":1,"zero_dm":0,"covariance":0,"phaseCovariance":0,"cal_source":0}
        self.subtractTipTilt=self.config.getVal("subtractTipTilt",default=0)
        self.phsphs=None
        self.storedPhs=None
        self.montePhaseCovMat=None
        self.subapVariance={}
        self.phssquared=None
        self.dmModeType=self.config.getVal("dmModeType")#specifies what the modes of the dm are - are they just pokes, or a more complicated shape fitted to the actuators?
        minarea=self.config.getVal("wfs_minarea")
        self.subflag=self.pupil.getSubapFlag(nsubx,minarea)#config.getVal("subflag")
        if self.dmflag!=None:#not a zonal DM.
            self.dmflag_1d=self.dmflag.ravel()#numpy.reshape(self.dmflag,(nact*nact,))
            self.dmindices=numpy.nonzero(self.dmflag_1d)[0]
            self.nacts=int(numpy.sum(self.dmflag))
        else:
            self.nacts=self.nact
        self.subflag_1d=self.subflag.ravel()#numpy.reshape(self.subflag,(nsubx*nsubx,))
        self.wfsdata = numpy.sum(self.subflag_1d)

        self.pxlscale=pxlToRadPistonTheoretical(nsubx,self.wfs_n,self.nimg,self.wfslambda,self.telDiam)#self.wfslambda/(self.telDiam/nsubx)/(self.nfft/self.wfs_n)#number of radians corresponding to centroid offset of 1 pixel.
        self.nthreads=self.config.getVal("nthreads",default="all")
        if self.nthreads=="all":
            self.nthreads=self.config.getVal("ncpu")
        self.monteNoiseCovariance=self.config.getVal("reconMonteCarloNoiseCovariance",default=1)
        self.noiseCov=self.config.getVal("reconNoiseCovariance",default=0.25)
        if self.dmModeType=="zernike":#use zernike modes but with a xinetics DM (ie actuators).
            #import util.tel
            #import util.zernikeMod
            pupiltmp=util.tel.Pupil(self.npup,self.npup*2,0)
            self.nmodes=self.config.getVal("dmModes",default=self.nacts)
            modes=numpy.array(util.zernikeMod.Zernike(pupiltmp,self.nmodes+1,computeInv=0).zern[1:,])/numpy.sqrt(numpy.sum(numpy.sum(self.pupil.fn)))
            util.zernikeMod.normalise(modes)
            fit=util.fitModes.fitModes(self.npup,self.nact,self.nmodes,modes,self.interpType,self.actCoupling,self.actFlattening)
            fit.fit()
            self.actModes=fit.actarr#this are used to fit actuators with control matrix.
            self.mirrorModes=fit.interpolatedAct.copy()#these are used to compute phase covariance.
            pfn=numpy.array(self.pupil.fn).astype(numpy.float32)
            for i in range(self.mirrorModes.shape[0]):
                self.mirrorModes[i]*=pfn
                #now give orthonormal scaling...
                self.mirrorModes[i]/=numpy.sqrt(numpy.sum(numpy.sum(self.mirrorModes[i]*self.mirrorModes[i])))
        elif self.dmModeType=="poke":
            self.nmodes=self.nacts
            self.mirrorModes=None
        elif self.dmModeType=="zernDM":#zernike modes on a zernike DM (science.zdm).  Here, the actuator values correspond to zernike values.
            self.mirrorModes=None
            self.nmodes=self.config.getVal("nmodes",default=int(self.nacts*0.75))
            self.nacts=self.nmodes
            self.Noll=numpy.zeros((self.nmodes,self.nmodes),numpy.float64)
            nollfile=string.join(__file__.split("/")[:-2]+["util","noll.fits"],"/")
            noll=util.FITS.Read(nollfile)[1]#the Noll matrix... doesn't include piston, so need to add.
            if noll.shape[0]+1<self.nmodes:
                    print "Noll matrix not large enough (noll.fits), please produce a larger one"
                    self.Noll[1::noll.shape[0]+1,1::noll.shape[1]+1]=noll
            else:
                    self.Noll[1:,1:]=noll[:self.nmodes-1,:self.nmodes-1]
            #correct the noll matrix for our diam and R0 (at 500nm) and then our wavelength....
            self.Noll*=(self.telDiam/self.r0)**(5./3)*(500./self.wfslambda)**2

            #convert the noll matrix to pixel scale.  Currently, units are radians of phase per radius of telescope, need to convert that into pixels in the centroided values.
            #pxlscale=self.wfslambda/(self.telDiam/nsubx)/(self.nfft/self.wfs_n)#number of radians corresponding to centroid offset of 1 pixel.  Possibly should be a factor of 1.22 larger?  But not scaled by 1e-9 since we want the result in nm anyway, and then arctan of this?
            #mphase=pxlscale*(self.telDiam/nsubx/2)#nano metres of phase intruduced by this tilt.
            #radpxl=mphase/self.wfslambda*2*numpy.pi#radians of phase which cause a shift of 1 pxl in centroid value.
            #radpxl*=nsubx#radians of phase across the entire pupil which give shift of 1 pxl in centroid value, ie radians per pxl.
            #This is now a scaling factor for the Noll matrix (radians^2 to pixels^2).
            #Now, dont scale to pixels, because the poke matrix scales from radians to pixels.
            self.Noll/=1#self.pxlscale**2#pxlToRadPistonTheoretical(nsubx,self.wfs_n,self.nfft,self.wfslambda,self.telDiam)  
        if forGUISetup:
            self.outputData=[(self.nacts,),numpy.float64]
        else:
            self.outputData=numpy.zeros((self.nacts,),numpy.float64)
            self.nsubaps=self.subflag.sum()
            self.fullWFSOutput=self.config.getVal("fullWFSOutput",default=1)
            if self.fullWFSOutput:
                self.centx=numpy.zeros((nsubx,nsubx),numpy.float64) #centroid arrays
                self.centy=numpy.zeros((nsubx,nsubx),numpy.float64)
            else:
                self.centx=numpy.zeros((self.nsubaps,),numpy.float64)
                self.centy=numpy.zeros((self.nsubaps,),numpy.float64)

	    self.ncents=numpy.sum(numpy.sum(self.subflag))*2
            self.pmxFilename=self.config.getVal("pmxFilename",default='poke_matrix.fits')#poke matrix will be saved as this.
            self.rmxFilename=self.config.getVal("rmxFilename",default=None,raiseerror=0)
            if self.reconType in ["regcg","regcg-waffle"]:#regularised conjugate gradient (Sophia di Moudi).
                self.regValue=self.config.getVal("regValue",default=0.02)#regularisation param for tichonov regularisation.
                self.regNIter=self.config.getVal("regNIter",raiseerror=0)#max number of iterations for cg for tichonov regularisation.
                self.lsreg=None
                self.regA=None
                self.regAFilename=self.config.getVal("regAFilename",raiseerror=0)#filename for regularised LHS.
                self.lsreg=util.regularisation.lsreg()

                if os.path.exists(self.pmxFilename):
                    self.pokemx=util.FITS.Read(self.pmxFilename)[1]
                    if self.pokemx.shape!=(self.nmodes,self.ncents):
                        raise Exception("pmx is wrong shape")
                else:
                    self.pokemx=numpy.zeros((self.nmodes,self.ncents),numpy.float64)
                if self.regAFilename!=None and os.path.exists(self.regAFilename):
                    self.regA=util.FITS.Read(self.regAFilename)[1]
                self.reconmx=None
            else:
                self.pokemx=numpy.zeros((self.nmodes,self.ncents),numpy.float32)
                self.reconmx=numpy.zeros((self.ncents,self.nmodes),numpy.float32)
            if self.rmxFilename!=None:
                if os.path.exists(self.rmxFilename):
                    if self.dmModeType=="zernDM":
                        self.reconmx[:,:]=util.FITS.Read(self.rmxFilename)[1]
                    else:
                        print "Loading reconmx from %s - shape expected %s"%(self.rmxFilename,str(self.reconmx.shape))
                        try:
                            self.reconmx[:]=util.FITS.Read(self.rmxFilename)[1]
                        except:
                            print "Unable to load %s - wrong shape?  Should be %s"%(self.rmxFilename,str(self.reconmx.shape))
                else:
                    print "WARNING - UNABLE TO FIND RECONSTRUCTOR %s"%self.rmxFilename
            self.tiltsens_list=[]#config.getVal("tiltsens_list")#this was tiltsens_list=['tiltsens4','tiltsens5','tiltsens6']
            self.cent_list=[]
            for k in self.parent.keys():#note, noise covariance stuff assumes only one cent parent.  As does other stuff, so this part is probably depreciated - ie only 1 cent parent is allowed...
                if k[:4]=="cent":
                    self.cent_list.append(k)
                elif k[:8]=="tiltsens":
                    self.tiltsens_list.append(k)
                elif k!="atmos":
                    print "Assuming parent %s is cent"%k
                    self.cent_list.append(k)
                    #raise Exception("xinterp_recon - unknown parent type %s"%k)
            self.poking=0
            # Set up count to only increase if subaperture is unvignetted (as defined by wfs_minarea)
            self.count=0
            self.lastControlCovariance=0
            self.Nph=self.config.getVal("wfs_sig")
            self.skybrightness=self.config.getVal("wfs_skybrightness")
            self.readNoiseSigma=self.config.getVal("wfs_read_sigma")
            self.readNoiseBg=self.config.getVal("wfs_read_mean")
            self.abortAfterPoke=self.config.getVal("abortAfterPoke",default=0)
            self.noiseMatrix=None
            self.phaseCov=None
            self.inputData=numpy.zeros(self.wfsdata*2,numpy.float64)
            self.decayFactor=self.config.getVal("decayFactor",default=1.)#a value used here when adding the new solution to the current actuator values.
            if self.reconType=="MAP":
                self.computeNoiseCovariance()
                self.computeControlMatrix=self.config.getVal("computeControlMatrix",default=0)
                if self.dmModeType=="zernDM":
                    self.phaseCov=self.Noll.copy()
                else:
                    if self.computeControlMatrix:
                        self.computePhaseCovariance()
                    else:
                        self.createMirrorModes()
            elif self.reconType=="pcg":
                #noisecov=self.computeSubapNoiseCovariance(self.wfs_n*self.wfs_n)
                if self.nact-1!=self.wfs_nsubx:
                    raise Exception("xinterp_recon - for PCG, nsubx must equal nact-1.")
                self.pcgCreateSparsePokeMx=self.config.getVal("pcgCreateSparsePokeMx",default=1)
                if self.pcgCreateSparsePokeMx==0:
                    self.pokemxPCGtmp=numpy.zeros((self.nact*self.nact,self.ncents),numpy.float64)
                    self.pokemxPCG=numpy.zeros((nsubx*nsubx*2,self.nact*self.nact),numpy.float64)
                else:
                    self.pokemxPCG=None

##                 if self.nacts!=self.nact*self.nact:
##                     print "WARNING: xinterp_recon - PCG won't work unless you have a full pokemx.  You need to set all of dmflag (in the config file) to 1"
##                     raise Exception("WARNING: xinterp_recon - PCG won't work unless you have a full pokemx.  You need to set all of dmflag (in the config file) to 1")
##                 if self.ncents!=self.nsubx*self.nsubx*2:
##                     print "WARNING: xinterp_recon - PCG won't work unless you have a full pokemx.  You need to set all of subflag (in the config file) to 1"
##                     raise Exception("WARNING: xinterp_recon - PCG won't work unless you have a full pokemx.  You need to set all of subflag (in the config file) to 1")
                if self.dmModeType!="poke":
                    print "WARNING: xinterp_recon - PCG will only work with a poked mirror (don't try to fit modes to it)"
                    raise Exception("WARNING: xinterp_recon - PCG will only work with a poked mirror (don't try to fit modes to it)")
                #self.pcg=util.fdpcg.createPcg(self.pokemxPCG,noisecov,self.nact,self.telDiam,self.r0,self.l0)
                print "xinterp_recon - creating centroid object for pixel factors"
                #import util.centroid#don't know why this is needed - python bug?
                c=util.centroid.centroid(self.wfs_nsubx,self.pupil.fn,readnoise=self.readNoiseSigma,addPoisson=1,noiseFloor=self.readNoiseSigma*3+self.readNoiseBg,readbg=self.readNoiseBg,sig=self.Nph)
                if self.monteNoiseCovariance:
                    self.noiseCov=c.computeNoiseCovariance(20,convertToRad=1)
                else:
                    pass
                    #self.noiseCov=#about right...
                self.radToPxlFactor=1/c.convFactor**2
                bccb=util.phaseCovariance.bccb(nact=self.nact,telDiam=self.telDiam,r0=self.r0,l0=self.l0,expand=0)
                print "xinterp_recon - Inverting the bccb phase covariance matrix (may take a while if large - except it won't because now using quick method) (shape=%s)"%str(bccb.shape)
                #self.phasecov=1./numpy.diagonal(inverse(bccb))
                ibccb=util.phaseCovariance.invertBCCB(bccb,self.nact)
                self.phasecov=numpy.zeros((bccb.shape[0],),numpy.float64)
                self.phasecov[:,]=1./ibccb[0]
                self.phasecov*=self.radToPxlFactor#to get into pixels...
                self.noiseCov*=self.radToPxlFactor#to get into pixels (from radians)
                self.pcgOptions=self.config.getVal("pcgOptions",default={})
                #self.atmosGeom=self.config.getVal("atmosGeom",default=None,raiseerror=0)
                if self.atmosGeom==None:#depreciated...
                    strLayer=self.config.getVal("strLayer")#depreciated
                    if type(strLayer)!=type([]):
                        strLayer=[strLayer]
                    nlayers=len(strLayer)
                else:
                    nlayers=len(self.atmosGeom.layerDict.keys())
                self.convergenceWarning=self.config.getVal("fdpcg_convergenceWarning",default=1)#if 1 will print a warning when convergence not reached.
                self.pcgGainFactor=self.config.getVal("pcgGainFactor",default=1.)#used to multiply the estimated phase by.
                self.pcgFakePokeMxVal=self.config.getVal("pcgFakePokeMxVal",default=0.31)#will use an artificial poke matrix if >0.  Probably best to do so.  The value here should be approximately the value that would be given by a real poke matrix, and can have a strong influence on strehl.
                if self.pcgFakePokeMxVal==0.31:
                    print "WARNING - you are using the default value of 0.31 for pcgFakePokeMxVal - this could mean your Strehl is less than it could be.  To get a better value, create a poke matrix and look at the typical values of poked centroids, and use this (and try varying a bit, eg maybe slightly less...)."
                if self.pcgFakePokeMxVal>0 and self.pcgCreateSparsePokeMx==0:
                    self.pokemxPCG[:,]=util.createPokeMx.straightPokes(self.nact,[self.pcgFakePokeMxVal]).transpose()
                
                self.pcgConvergenceValue=self.config.getVal("pcgConvergenceValue",default=0.1)#set this lower (eg 0.001) for slight improvement in reconstruction.  This value is the maximum allowed difference between the pcg solution from the previous iteration, before it is considered to have converged.
                self.pcgMaxIter=self.config.getVal("pcgMaxIter",default=50)
                print "xinterp_recon - creating pcg object"
                self.pcg=util.fdpcg.pcg(self.pokemxPCG,self.noiseCov,self.phasecov,None,nact=self.nact,turbStrength=range(nlayers),options=self.pcgOptions,minIter=1,maxIter=self.pcgMaxIter,convergenceValue=self.pcgConvergenceValue,convergenceWarning=self.convergenceWarning,gainfactor=self.pcgGainFactor,fakePokeMxVal=self.pcgFakePokeMxVal)

    def newParent(self,parent,idstr=None):
        if type(parent)!=type({}):
            parent={"cent":parent}
        self.parent=parent
        self.tiltsens_list=[]#config.getVal("tiltsens_list")#this was tiltsens_list=['tiltsens4','tiltsens5','tiltsens6']
        self.cent_list=[]
        for k in self.parent.keys():
            if k[:4]=="cent":
                self.cent_list.append(k)
            elif k[:8]=="tiltsens":
                self.tiltsens_list.append(k)
            elif k!="atmos":
                print "Assuming parent %s is cent"%k
                self.cent_list.append(k)
                #raise Exception("xinterp_recon - unknown parent type %s"%k)
                
### Recontructor main loop #############################################


    def generateNext(self,msg=None):
        """Data coming in from parents can be of 2 types:
         - centroid data (lgs)
         - tilt sensitivity data (ngs)
        Determine which is which from the name of the parent object
        (the dictionary key).  If cent*, is centroid data, if
        tiltsens*, is tilt sensitivity data.
        """
        if self.generate==1:
            if self.newDataWaiting:
                nin=0
                for key in self.parent.keys():
                    if self.parent[key].dataValid==1:
                        nin+=1
                    else:
                        print "Xinterp Recon: Waiting for data from wfs, but not valid"
                if nin>0:
                    if nin==len(self.parent.keys()):
                        self.dataValid=1
                    else:
                        print "Xinterp Recon: Warning - got some data but not all, setting dataValid=0"
                        self.dataValid=0
                else:
                    self.dataValid=0
            if self.dataValid:
                if self.control["covariance"]:
                    self.noiseCovariance()
                else:
                    self.lastControlCovariance=0
                if self.parent.has_key("atmos") and self.control["phaseCovariance"]:
                    self.montePhaseCovariance()
                    self.computeMeanPhaseVariance()
                self.calc()# Calculate DM updates
        else:
            self.dataValid=0
    def calc(self):
        nsubx=self.wfs_nsubx # Number of subapertures across pupil
        wfsdata=self.wfsdata # Number of flagged subapertures
        pokemx=self.pokemx
        if self.control["poke"]:
            self.poking=1
            self.control['poke']=0
            self.pokemx[:,]=0
            self.control["close_dm"]=0
            print "Will be poking for %d iterations"%(self.nmodes+1)
        if self.control["zero_dm"]:
            self.control["zero_dm"]=0
            self.outputData[:,]=0.
        self.centx*=0.#numpy.zeros((nsubx,nsubx),numpy.float64)# Centroid arrays
        self.centy*=0.#numpy.zeros((nsubx,nsubx),numpy.float64)

	if len(self.tiltsens_list) > 0:
	    self.tiptilt*=0.#numpy.zeros(2,numpy.float64)# Closed loop tilt increments
            for tiltsensproc in self.tiltsens_list:
                self.tiptilt+=self.parent[tiltsensproc].outputData*4.0#new centroids
                # 4 is for scaling from fsm module to Zernike tip/tilt
            self.tiptilt /= float(len(self.tiltsens_list))
        for centproc in self.cent_list:  #centproc_list:
            if self.fullWFSOutput:
                self.centx += self.parent[centproc].outputData[:,:,0]#tempcentx
                self.centy += self.parent[centproc].outputData[:,:,1]#tempcenty
            else:
                self.centx+=self.parent[centproc].outputData[:,0]
                self.centy+=self.parent[centproc].outputData[:,1]
        #Now average the centroids
        self.centx /= float(len(self.cent_list))  #centproc_list))
        self.centy /= float(len(self.cent_list))  #centproc_list))
        if self.subtractTipTilt==-1 or (self.subtractTipTilt==1 and self.control["cal_source"]==0 and self.control["poke"]==0 and self.poking==0):
            N=self.nsubaps
            self.centx-=self.centx.sum()/N
            self.centy-=self.centy.sum()/N
            
        #print self.centx.shape,self.centy.shape,self.subflag_1d.shape
        # self.calc2()
        if self.removeAverageCent:#probably not the correct way of doing it since centx is all centroids, not just the used ones.
            self.centx-=numpy.average(self.centx)
            self.centy-=numpy.average(self.centy)
        if self.poking>0 and self.poking<=self.nmodes:
            #First, set the next actuator to be poked... note, this only bothers with actuators that have some effect...
            if self.dmModeType=="poke":
                self.outputData[:,]=0.
                if self.reconType=="MAP":
                    pokeval=self.pokeval/self.mirrorScale[self.poking-1]
                else:
                    pokeval=self.pokeval
                self.outputData[self.poking-1]=pokeval
            elif self.dmModeType=="zernike":
                self.outputData[:,]=numpy.take(self.actModes[self.poking-1].ravel(),self.dmindices)*self.pokeval#place the mode onto the actuators.
            elif self.dmModeType=="zernDM":
                self.outputData[:,]=0.
                self.outputData[self.poking-1]=self.pokeval
        if self.poking>1 and self.poking<=self.nmodes+1:
            #then use the centroid values from previous poked actuator
            #to fill the poke matrix.
            if self.fullWFSOutput:
                pokemx[self.poking-2,:self.wfsdata]=numpy.compress(self.subflag_1d!=0,numpy.reshape(self.centx,(nsubx*nsubx,)))#Uses only flagged subapertures
                pokemx[self.poking-2,self.wfsdata:,]=numpy.compress(self.subflag_1d!=0,numpy.reshape(self.centy,(nsubx*nsubx,)))#Uses only flagged subapertures
            else:
                pokemx[self.poking-2,:self.wfsdata]=self.centx
                pokemx[self.poking-2,self.wfsdata:,]=self.centy
                
        if self.poking==self.nmodes+1:#flatten the dm
            self.outputData[:,]=0.
        if self.poking==self.nmodes+1:#finished poking...
            print "Calculating control matrix"
            #pokemx[:,] = numpy.where(pokemx>0.0,pokemx,0) #Ignore values that are less than 0.01
            start=0
            pokemx/=self.pokeval#*self.gamma # scale pokematrix by poke value.  Gain value removed on 080723 by agb.
            util.FITS.Write(pokemx,self.pmxFilename)
            if self.dmModeType=="zernDM":
                pokemx=pokemx[1:]#remove piston.
                start=1
                if self.reconmx!=None:
                    self.reconmx*=0
            if self.reconType=="leastSquares" or self.reconType=="pinv":
                self.reconmx[:,start:]=numpy.linalg.pinv(pokemx,self.rcond)
            elif self.reconType=="regBig":#big inversion
                self.reconmx[:,start:]=util.regularisation.invert(pokemx,self.rcond,large=1)
            elif self.reconType=="regSmall" or self.reconType=="reg" or self.reconType=="regularised":#small inversion
                self.reconmx[:,start:]=util.regularisation.invert(pokemx,self.rcond,large=0)
            elif self.reconType=="SVD":
                self.reconmx[:,start:]=self.control_matrix(pokemx,self.rcond)
            elif self.reconType=="MAP":
                #if self.dmModeType=="poke":#scale the modes by mirrorScale (since this was also used for poking).
                #    for i in range(self.nmodes):
                #        pokemx[i]*=self.mirrorScale[i]
                if self.computeControlMatrix:
                    self.reconmx[:,start:]=self.map_matrix(pokemx,self.phaseCov,self.noiseMatrix,self.rcond)
            elif self.reconType=="pcg":
                #pcg is not a MVM operator...
                #print "TODO: beef pokemx up to include all modes and cents"
                pokemx/=self.gamma # scale pokematrix by poke value
                self.expandPmx()
                if self.pcgFakePokeMxVal>0 and self.pcgCreateSparsePokeMx==0:
                    self.pokemxPCGOld=self.pokemxPCG.copy()
                    self.pokemxPCG[:,]=util.createPokeMx.straightPokes(self.nact,[self.pcgFakePokeMxVal]).transpose()#,self.pupil.subflag,returnfull=1))
                self.pcg.newPokemx(self.pokemxPCG)
            elif self.reconType=="regcg":#regularised conjugate gradient (Sophia di Moudi).
                self.regA=self.lsreg.regmatrix1(pokemx,self.regValue)
            elif self.reconType=="regcg-waffle":#regularised least squares with waffle penalisation cg.
                self.regA=self.lsreg.wafflemx(pokemx,self.regValue)
            if self.reconmx!=None:
                util.FITS.Write(self.reconmx,self.rmxFilename)
            self.poking = 0 #Finish poking
            if self.abortAfterPoke:
                print "Finished poking - aborting simulation"
                Scientific.MPI.world.abort(0)#this raises an error if python, or aborts correctly if mpipython - however, result is as desired - the simulation finishes.
        if self.poking>0:
            self.poking+=1 #Move to next actuator
        if self.control["close_dm"]:#do the reconstruction.
            self.calc2()
        #Now forward stuff to the DM.
        #self.outputData=self.coeff2
    def expandPmx(self):
        """expand self.pokemx into self.pokemxPCG which includes all centroids and actuators."""
        actpos=0
        for i in range(self.nact):
            for j in range(self.nact):
                if self.dmflag[i,j]==1:
                    self.pokemxPCGtmp[i*self.nact+j]=self.pokemx[actpos]
                    actpos+=1
        cpos=0
        nsa=self.wfs_nsubx*self.wfs_nsubx
        nc=self.ncents/2
        for i in range(self.wfs_nsubx):
            for j in range(self.wfs_nsubx):
                if self.subflag[i,j]==1:
                    self.pokemxPCG[i*self.wfs_nsubx+j,:,]=self.pokemxPCGtmp[:,cpos]
                    self.pokemxPCG[i*self.wfs_nsubx+j+nsa,:,]=self.pokemxPCGtmp[:,cpos+nc]
                    cpos+=1
        
    def calc2(self):
        """Reconstruct using poke matrix or pcg or CG (Added regularised CG version)"""
        nsubx=self.wfs_nsubx
        wfsdata=self.wfsdata
        data=self.inputData#numpy.zeros(wfsdata*2,numpy.float64)
        if self.fullWFSOutput:
            data[:wfsdata]=numpy.compress(self.subflag_1d!=0,numpy.reshape(self.centx,(nsubx*nsubx,)))#Uses only flagged subapertures
            data[wfsdata:,]=numpy.compress(self.subflag_1d!=0,numpy.reshape(self.centy,(nsubx*nsubx,)))#Uses only flagged subapertures
        else:
            data[:wfsdata]=self.centx
            data[wfsdata:]=self.centy
        #self.outputData*=(1-self.gamma)#this line makes performance worse.
        if self.reconType=="pcg":
            self.pcg.solve(data,usePrevious=0)#for some reason, usePrevious=1 will cause it to blow up eventually... have no idea why - maybe the pcg algorithm is not too good at starting points close to the initial.
            #print "TODO: select only needed phase values - only the used acts"
            self.outputData[:,]+=-self.pcg.gainfactor*numpy.take(numpy.array(self.pcg.x),self.dmindices)
        elif self.reconType in ["regcg","regcg-waffle"]:
            self.regb=numpy.dot(self.pokemx,data)#dot poke matrix with slopes to get RHS.
            self.outputData*=self.decayFactor
            if self.regA==None:
                print "Computing regularised A matrix"
                self.regA=self.lsreg.regmatrix1(self.pokemx,self.regValue)
            self.outputData[:]+=-self.gamma*self.lsreg.solvecg(self.regA,self.regb,maxiter=self.regNIter)
        else:
            if self.dmModeType=="poke":
                self.outputData*=self.decayFactor
                self.outputData[:,]+=-self.gamma*numpy.dot(data,self.reconmx)
            elif self.dmModeType=="zernike":
                modes=numpy.dot(data,self.reconmx)
                output=numpy.zeros((self.nact,self.nact),numpy.float64)
                for i in range(self.nmodes):
                    output+=self.actModes[i]*modes[i]
                self.outputData[:,]+=-self.gamma*numpy.take(output.ravel(),self.dmindices)
            elif self.dmModeType=="zernDM":
                self.outputData*=self.decayFactor
                self.outputData+=-self.gamma*numpy.dot(data,self.reconmx)

    def control_matrix(self,pokemx,corr_thresh):
        """Make control matrix from poke matrix"""
        svec=numpy.linalg.svd(pokemx)
        u,a,vt=svec

        util.FITS.Write(a,'singular_values.fits')
        #a[self.n_modes_corrected:,]=0.
        #print "WARNING WARNING WARNING: xinterp_recon doing SVD control matrix currently zeros elements if they are greater than corr_thresh.  Shouldn't this be if they are less than this???  Please investigate and fix sometime..."
        n_removed=0
        for i in range(len(a)):
            if a[i] < corr_thresh:
                print "removing eval %d (value %g)"%(i,a[i])
                a[i]=0.
                n_removed += 1
        print 'Removed %d modes from control matrix (threshold %g)'%(n_removed,corr_thresh)
        print 'Eigenvalues:',a

        ut=numpy.transpose(u)
        v=numpy.transpose(vt)
        id=numpy.identity(a.shape[0])
        ai=numpy.multiply(a,id)
        ai=numpy.where(ai != 0, 1/ai, 0)
        #print ai.shape,ut.shape,pokemx.shape
        tmp=numpy.dot(ai[:,:ut.shape[0]], ut[:ai.shape[1]])
        reconmx = numpy.dot(v[:,:tmp.shape[0]], tmp)
        return reconmx

    def map_matrix(self,pokemx,phaseCov,noiseCov,rcond):
        """Make a control matrix using MAP information"""
        #note, phaseCov would be the noll matrix (scaled) if zernike modes.
        #noise matrix is the subaperture noise covariance matrix.
        mapVersion=0#change this to test different map versions
        pmx=numpy.transpose(pokemx)
        pmxT=pokemx
        if mapVersion==0:
            # This was the original version - not sure if correct...
            reconmx=numpy.transpose(numpy.dot(numpy.dot(phaseCov,pmxT),numpy.linalg.pinv(numpy.dot(numpy.dot(pmx,phaseCov),pmxT)+noiseCov,rcond)))
        elif mapVersion==1:
            # Actually, maybe this should be (Roggerman, imaging through turbulence book, pg 193):
            invNoiseCov=noiseCov.copy()
            diag=noiseCov[::noiseCov.shape[0]+1]
            invNoiseCov[::noiseCov.shape[0]+1]=numpy.where(diag==0,0,1./diag)#diagonal anyway
            invPhaseCov=numpy.linalg.inv(phaseCov)
            reconmx=numpy.dot(numpy.linalg.inv(numpy.dot(pmxT,numpy.dot(invNoiseCov,pmx))+invPhaseCov),numpy.dot(pmxT,invNoiseCov))
            reconmx=numpy.transpose(reconmx)#transpose since we do the dot in a funny order
        if self.dmModeType=="poke":
            #need to scale the xinetics_dm output to unity...
            #(for MAP, we are assuming that mirror modes (pokes) are scaled s.t. int(mode**2)==1.
            #we can do this by altering the reconmx here.
            for i in range(self.nmodes):
                reconmx[:,i]/=self.mirrorScale[i]

            
        return reconmx
    
    def createMirrorModes(self,actmodeArr=None,fitpup=1):
        """Fit the mirror to the actuators."""
        global util#why is this needed?  python bug?
        if self.dmModeType=="poke":
            self.mirrorScale=numpy.zeros((self.nmodes),numpy.float64)
            createActs=0
            self.mirrorModes=numpy.zeros((self.nmodes,self.npup,self.npup),numpy.float32)#this are the "modes", which might be individual actuators, or might be eg zernikes forced on the mirror etc.
            actmode=numpy.zeros((self.nact,self.nact),numpy.float32)
            if type(actmodeArr)==type(None):
                createActs=1
            else:
                actmodeArr=actmodeArr.astype(numpy.float32)
            pfn=numpy.array(self.pupil.fn).astype(numpy.float32)
            print "xinterp_recon: fitting mirror modes"
            mirrorSurface=MirrorSurface(self.interpType,self.npup,self.nact,actoffset=self.actoffset,actCoupling=self.actCoupling,actFlattening=self.actFlattening)
            for i in range(self.nmodes):
                y=self.dmindices[i]/self.nact
                x=self.dmindices[i]%self.nact
                if createActs:
                    actmode[y,x]=1.
                else:
                    numpy.put(actmode.ravel(),self.dmindices,actmodeArr[i])

                
                #This assumes that we are using a xinterp_dm module...
                self.mirrorModes[i]=mirrorSurface.fit(actmode)#util.fitModes.interp(actmode,self.interpType,self.npup,self.actCoupling,self.actFlattening,self.nact)
                if fitpup:
                    self.mirrorModes[i]*=pfn

                #now give orthonormal scaling... (even if not ortho basis).
                self.mirrorScale[i]=numpy.sqrt(numpy.sum(self.mirrorModes[i]*self.mirrorModes[i]))
                self.mirrorModes[i]/=self.mirrorScale[i]
                #and then apply the pupil function.  Not sure this is desirable...
                #self.mirrorModes[i]*=pfn
                if createActs:
                    actmode[y,x]=0.
        elif self.dmModeType=="zernDM":
            print "create mirror modes for zernDM: using zernikes..."
            global util
            self.mirrorModes=util.zernikeMod.Zernike(self.npup,self.nmodes,computeInv=0).zern.astype(numpy.float32)#keep piston on purpose - it is ignored later (ie the zdm ignores it).  Kept for historical reasons, and probably for completeness.
            util.zernikeMod.normalise(self.mirrorModes)
        else:
            raise Exception("Mirror modes already created - xinterp_recon")

    def makeMirrorModes(self):
        """Takes the poke matrix, performs SVD, and then uses the mirror
        mode matrix created to make the fundamental mirror modes.
        This is typically called by the GUI, which also called the method in
        the DM object too, after passing it the modeMatrix."""
        u,a,vt=numpy.linalg.svd(self.pokemx)
        return u#this can be used with createMirrorModes(makeMirrorModes()) to create the natural mirror modes... not sure how useful that is though!
    def computeNoiseCovariance(self):
        """Compute the centroid/subap noise covariance.  This was taken
        directly from the map reconstructor in glao_zernike_recon."""
        print "Computing noise covariance matrix"
        self.noiseMatrix=numpy.zeros((self.ncents,self.ncents),numpy.float64)#noise matrix (diagonal) in radians of phase.
        nsubx=self.wfs_nsubx
        wfsn=self.wfs_n
        i=0
        for y in range(nsubx):
            for x in range(nsubx):
                if self.subflag[y,x]:
                    # N. of photons in this subap (edge of pupil?)  should it include wfs_ncen ?  Probably...
                    Npxl=numpy.sum(numpy.sum(self.pupil.fn[y*wfsn:y*wfsn+wfsn,x*wfsn:x*wfsn+wfsn]))
                    ncv=self.computeSubapNoiseCovariance(Npxl)
                    self.noiseMatrix[i,i]=ncv
                    self.noiseMatrix[i+self.ncents/2,i+self.ncents/2]=ncv
                    i+=1
                    #Note, noiseMatrix is in radians^2 of phase.  So convert into pixels.
        self.noiseMatrixtheta=self.noiseMatrix.copy()#save for reference...
        #theta=self.noiseMatrix/numpy.pi**2*self.pxlscale#FA suggests this...
        #self.noiseMatrix=theta/self.pxlscale#ie basically divide by pi^2.
        self.noiseMatrix/=self.pxlscale**2
        self.noiseMatrixnew=self.noiseMatrix.copy()
    def computeSubapNoiseCovariance(self,Npxl):
        """Compute noise covariance of centroid measurements.
        Note - this could well be wrong..."""
        wfsn=self.wfs_n
        Nt=1.*(self.nfft/self.wfs_n)
        Nd=Nt#Note, Nt is the FWHM, and Nd is FWHM of diffraction pattern.  Nt will be bigger due to atmospheric smearing.  
        
        Nph=self.Nph*Npxl/wfsn**2
        # N of photons in sky background
        Nbg=self.skybrightness*Npxl
        # Compute the SH measurement error due to photon noise (eq 5.39 in Rodier book)
        if Nph==0:
            Nph=1e-6
        self.SHmeasErrPhoton=numpy.pi**2/2./Nph*(Nt/Nd)**2
        # and due to readout noise: eq 5.40 in Rodier book
        #self.SHmeasErrRead=((numpy.pi*self.readNoiseSigma*Npxl/Nph/Nd)**2)/3.#This doesn't seem to give correct results.
        #Try agb version, for better
        #agreement with monte-carlo...
        #note, this has no theoretical reason
        #for being.
        self.SHmeasErrRead=self.readNoiseSigma/Nph/Nd*numpy.pi/3*10
        #and due to sky background: eq 5.41 in Rodier
        self.SHmeasErrBg=4/3.*(numpy.pi/Nph*Nt/Nd)**2*Nbg
        return self.SHmeasErrPhoton+self.SHmeasErrRead+self.SHmeasErrBg

    def computePhaseCovariance(self):
        """Compute phase covariance of the mirror modes"""
        #first check to see if one has been saved previously...
        myfilename="phaseCov_%s_%d_%d_%g_%g_%g_%g"%("vk",self.npup,self.nmodes,self.r0,self.telDiam,self.l0,self.dm.reconLam)
        lenmyfilename=len(myfilename)
        dirlist=os.listdir("/tmp")
        dataToUse=None
        for filename in dirlist:
            if filename[:lenmyfilename]==myfilename:
                print "possible phase covariance file found: %s"%filename
                data=util.FITS.Read("/tmp/"+filename)
                if len(data)==6:#head, covariance, head, pupil, head, modes
                    if numpy.sum(data[3]==self.pupil.fn)==self.pupil.fn.shape[0]*self.pupil.fn.shape[1]:
                        #pupils agree
                        if type(self.mirrorModes)==type(None):
                            self.createMirrorModes()#make the modes...
                        if numpy.sum(data[5]==self.mirrorModes)==self.mirrorModes.shape[0]*self.mirrorModes.shape[1]*self.mirrorModes.shape[2]:
                            #modes also agree...
                            if data[0]["parsed"].has_key("COVERSN"):
                                if data[0]["parsed"]["COVERSN"]==util.phaseCovariance.version():
                                    #version numbers agree...
                                    # so we can use this file.
                                    dataToUse=data[1]
                                else:
                                    print "phase covariance file has wrong version"
                            else:
                                print "phase covariance file has no version"
                        else:
                            print "phase covariance file has wrong modes"
                    else:
                        print "phase covariance file has wrong pupil"
                else:
                    print "phase covariance file contains wrong FITS extensions"
            if type(dataToUse)!=type(None):
                break
        if type(dataToUse)==type(None):
            #generate the matrix...
            #o=self.dm.computePhaseCovariance(self.atmosGeom,self.pupil.r2,self.r0,self.l0,self.dm.reconLam,nthreads=self.nthreads,mirrorSurface=self.mirrorSurface,width=self.mirrorModeWidth)
            if type(self.mirrorModes)==type(None):
                self.createMirrorModes()#make the modes
            p,m,v,o=util.phaseCovariance.make(typ="vk",npup=self.npup,nmode=self.nmodes,pupil=self.pupil,modes=self.mirrorModes,r0=self.r0,telDiam=self.telDiam,l0=self.l0,nthreads=self.nthreads)
            self.phaseCov=o
            #and then save if...
            num=""
            cnt=0
            while os.path.exists("/tmp/"+myfilename+num+".fits"):
                cnt+=1
                num="%03d"%cnt
            myfilename="/tmp/"+myfilename+num+".fits"
            print "Saving phase covariance data to file %s"%myfilename
            util.FITS.Write(self.phaseCov,myfilename,extraHeader=["COVERSN = %s"%str(util.phaseCovariance.version())])
            util.FITS.Write(self.pupil.fn,myfilename,writeMode="a")
            util.FITS.Write(self.mirrorModes,myfilename,writeMode="a")
        else:
            print "Using phase covariance data from file %s"%filename
            self.phaseCov=dataToUse
        # this is the phase covariance matrix, probably in radians.
        # correct the phase matrix for our diam and R0.
        #self.phaseCov*=(self.telDiam/self.r0)**(5./3)
        #convert the noll matrix to pixel scale.  Currently, units are radians of phase per radius of telescope, need to convert that into pixels in the centroided values.
        #Note, self.pxlscale=self.wfslambda/(self.telDiam/nsubx)/(self.nfft/self.wfs_n)#number of radians corresponding to centroid offset of 1 pixel
        if self.dmModeType=="poke":
            #scaling of phaseCov is already correct.
            #could try multiplying by gamma...though this is probably incorrect
            #self.phaseCov*=self.gamma*self.gamma
            pass
        elif self.dmModeType=="zernike":
            #mphase=self.pxlscale*(self.telDiam/self.wfs_nsubx/2)#metres of phase intruduced by this tilt.
            #radpxl=mphase/self.wfslambda*2*numpy.pi#radians of phase which cause a shift of 1 pxl in centroid value.
            #radpxl*=self.wfs_nsubx#radians of phase across the entire pupil which give shift of 1 pxl in centroid value, ie radians per pxl.
            # radpxl is pi * wfs_n * wfs_nsubx/nfft
            # This is now a scaling factor for the phaseCov matrix.
            self.phaseCov/=1#self.pxlscale**2
        elif self.dmModeType=="zernDM":
            #print "TODO: dmModeType==zernDM, computePhaseCovariance"
            self.phaseCov/=1#self.pxlscale**2#convert from radians^2 to pixels^2.
        #self.phaseCov=numpy.dot(self.invModeInteraction,self.phaseCov)
    def noiseCovariance(self):
        """This can be used to compute the wfs noise covariance matrix
        - if the DM output is zero... ie noiseless..."""
        if self.lastControlCovariance==0:#was zero last time - so start fresh.
            self.sumcent={}
            self.sum2cent={}
            self.subapVariance={}
            self.covcnt=0
            #for key in self.parent.keys():
            key="cent"
            self.sumcent[key]=numpy.zeros(self.parent[key].outputData.shape,self.parent[key].outputData.dtype)
            self.sum2cent[key]=numpy.zeros(self.parent[key].outputData.shape,self.parent[key].outputData.dtype)
        else:
            self.covcnt+=1
            #for key in self.parent.keys():
            key="cent"
            self.sumcent[key]+=self.parent[key].outputData
            self.sum2cent[key]+=self.parent[key].outputData*self.parent[key].outputData
        #for key in self.parent.keys():
        key="cent"
        self.subapVariance[key]=self.sum2cent[key]/self.covcnt-self.sumcent[key]*self.sumcent[key]/(self.covcnt*self.covcnt)
        self.lastControlCovariance=self.control["covariance"]


    def getNoiseCovariance(self,name="centroidNoiseCovariance"):
        """This can be called when a science.centCovariance is also being used and has computed covariances.
        Typically, will be called from the GUI.
        """
        self.subapVariance["cent"]=self.config.postGet(name)

    def makeModeInteractionMatrix(self):
        n=self.mirrorModes.shape[0]
        imat=numpy.zeros((n,n),numpy.float64)
        for i in range(n):
            for j in range(i,n):
                imat[i,j]=numpy.sum(self.mirrorModes[i]*self.mirrorModes[j])
                imat[j,i]=imat[i,j]
        self.modeInteractionMatrix=imat
        return imat

        
    def makeMonteRecon(self,mapVersion=0):
        """This can be called by GUI to make a montecarlo reconstructor
        using the montecarlo noise covariance.
        Makes a MAP reconstructor using phase and noise covariance matricees.
        """
        nsubap=self.wfs_nsubx*self.wfs_nsubx
        i=0
        if self.subapVariance.has_key("cent"):
            print "Using monte-noise covariance"
            for y in range(self.wfs_nsubx):
                for x in range(self.wfs_nsubx):
                    if self.subflag[y,x]:
                        self.noiseMatrix[i,i]=self.subapVariance["cent"][y,x,0]
                        self.noiseMatrix[i+self.ncents/2,i+self.ncents/2]=self.subapVariance["cent"][y,x,1]
                        i+=1
        if type(self.montePhaseCovMat)==type(None):
            print "Multiplying phasecov by pupil.sum"
            phaseCov=self.phaseCov*self.pupil.sum#this scaling factor is needed for zernDM since the zernikes are orthonormal, not in units of radians of phase.
        else:
            print "Using monte-phase covariance"
            #print "Multiplying phasecov by pupil.sum"
            #phaseCov=self.montePhaseCovMat*self.pupil.sum
        pmx=numpy.transpose(self.pokemx)
        pmxT=self.pokemx
        if mapVersion==0:
            self.reconmx=numpy.transpose(numpy.dot(numpy.dot(phaseCov,pmxT),numpy.linalg.pinv(numpy.dot(numpy.dot(pmx,phaseCov),pmxT)+self.noiseMatrix,self.rcond)))
            #Note, this is okay for matrix inversion lemma, provided noiseMatrix can be inverted easily (ie if diagonal or sparse).  Then, the inverse is of size nacts,nacts, rather than ncents,ncents.


        elif mapVersion==1:
            # Actually, maybe this should be (Roggerman, imaging through turbulence book, pg 193):  Indeed, these should be the same - and do seem to be - this could be a test of the relative scaling - to see whether it is correct! (remember to take mirrorScale into account if doing this in an interactive session).
            #The inversion here is size nacts,nacts, as is the invPhaseCov.  The invNoiseCov is size ncents,ncents, but usually diagonal, so okay.
            invNoiseCov=self.noiseMatrix.copy()
            diag=noiseCov[::noiseCov.shape[0]+1]
            invNoiseCov[::noiseCov.shape[0]+1]=numpy.where(diag==0,0,1./diag)#diagonal anyway
            invPhaseCov=numpy.linalg.inv(phaseCov)
            reconmx=numpy.dot(numpy.linalg.inv(numpy.dot(pmxT,numpy.dot(invNoiseCov,pmx))+invPhaseCov),numpy.dot(pmxT,invNoiseCov))
            self.reconmx=numpy.transpose(reconmx)#transpose since we do the dot in a funny order

        if self.dmModeType=="poke":
            #need to scale the xinetics_dm output to unity...
            #(for MAP, we are assuming that mirror modes (pokes) are scaled s.t. int(mode**2)==1.
            #we can do this by altering the reconmx here.
            for i in range(self.nmodes):
                self.reconmx[:,i]/=self.mirrorScale[i]


        print "Computed reconstructor (if monte-carlo noise and phase covariances have previously been computed, these will have been used), saved as monteReconmx.fits"
        util.FITS.Write(self.reconmx,"monteReconmx.fits")
        util.FITS.Write(phaseCov,"monteReconmx.fits",writeMode="a")
        util.FITS.Write(self.noiseMatrix,"monteReconmx.fits",writeMode="a")
        util.FITS.Write(pmxT,"monteReconmx.fits",writeMode="a")
        util.FITS.Write(self.mirrorScale,"monteReconmx.fits",writeMode="a")
    def createMirrorModeArr(self):
        """Here, creates the mirror modes, which may be actuators forced to zernike, or may be just the actuators themselves."""

    def constructPhase(self,phs=None):
        """Decompose the phase into the current mirror modes, and then reconstruct using these.
        Typically, this may be used by the GUI for inspection purposes... to see how orthogonal a basis is etc, ie how well phase can be represented in a given set of modes.
        """
        if type(phs)==type(None):
            phs=self.parent["atmos"].outputData
        a=numpy.zeros((self.nmodes,),numpy.float64)
        output=numpy.zeros(phs.shape,phs.dtype)
        if type(self.mirrorModes)==type(None):
            self.createMirrorModes()
        start=0
        if self.dmModeType=="zernDM":
            start=1
        for i in xrange(start,self.nmodes):
            a[i]=numpy.sum(phs*self.mirrorModes[i])
        for i in xrange(start,self.nmodes):
            output+=a[i]*self.mirrorModes[i]
        return output
    
    def montePhaseCovariance(self):
        a=numpy.zeros((self.nmodes,),numpy.float64)
        phs=self.parent["atmos"].outputData
        if type(self.mirrorModes)==type(None):
            self.createMirrorModes()
        start=0
        if self.dmModeType=="zernDM":
            start=1
        for i in xrange(start,self.nmodes):
            a[i]=numpy.sum(phs*self.mirrorModes[i])
        if type(self.monteNoll)==type(None):
            self.monteNoll=numpy.zeros((self.nmodes,self.nmodes),numpy.float64)
            self.navNoll=0
        self.monteNoll+=a[:,numpy.newaxis]*a[numpy.newaxis,:]
        self.navNoll+=1
    def finishMontePhaseCovariance(self,fname="montePhaseCovMat.fits"):
        self.montePhaseCovMat=self.monteNoll/(self.navNoll)#*self.pupil.sum)#*self.pxlscale**2)#should pupil.sum be with or without the central obscuration? (not that it makes a huge difference!).  Don't scale to pixels, because the poke matrix scales from radians to pixels.
        if fname!=None:
            print "%s has been written"%fname
            util.FITS.Write(self.montePhaseCovMat,fname)
        #This should be similar to self.Noll (check scaling - am not sure its quite right).
    def montePhaseCovarianceOrig(self):
        """Start computing the phase covariance from monte carlo.  This assumes have access to the atmospheric phase (unphysical), and adds the phase for many iterations, ready to be used by finishMontePhaseCovariance.
        """
        phs=self.parent["atmos"].outputData
        #First, want to create the <phi phi> array.  Since this is basically just the difference between two phases, can do that directly from this array.
        if type(self.storedPhs)==type(None):
            self.storedPhs=phs*phs[0,0]
            self.nPhs=1
        else:
            self.storedPhs+=phs*phs[0,0]
            self.nPhs+=1
        #result can then be used by finishMontePhaseCovariance().
        #note, this is actually the phase covariance, but what we're after is the phase covariance wrt the mirror modes.
        self.phsphs=self.storedPhs/self.nPhs#this is then the variance I think, since the mean should be zero.
        
    def finishMontePhaseCovarianceOrig(self):
        """Here, we compute the phase covariance from monte carlo..."""
        covMat=numpy.zeros((self.nmodes,self.nmodes),numpy.float32)
        if type(self.mirrorModes)==type(None):
            self.createMirrorModes()
        #not sure this is the correct thing to do:???
        cmod.phaseCov.covariance(self.phsphs.astype(numpy.float32),self.mirrorModes,covMat)
        self.montePhaseCovMat=covMat/(numpy.sum(self.pupil.fn))#*self.pxlscale**2)#convert from radians^2 to pixels^2.
        self.montePhaseCovMat[0]=0.#remove piston
        self.montePhaseCovMat[:,0]=0.#remove piston
        
        util.FITS.Write(self.montePhaseCovMat,"montePhaseCovMat.fits")
        return covMat
    def computeMeanPhaseVariance(self):
        """Compute the variance of phase at an arbitrary point..."""
        phs=self.parent["atmos"].outputData
        if type(self.phssquared)==type(None):
            self.phssquared=phs*phs
            self.phssummed=phs.copy()
            self.nPhsVar=1
        else:
            self.phssquared+=phs*phs
            self.phssummed+=phs
            self.nPhsVar+=1
        mean=self.phssummed/self.nPhsVar
        self.phaseVariance=self.phssquared/self.nPhsVar-mean*mean
        self.avPhaseVariance=numpy.sum(self.phaseVariance)/(self.npup*self.npup)


    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr==None or self.idstr==[None]:
            id=""
        else:
            id=" (%s)"%self.idstr
        txt=""
        if self.reconType=="pcg":
            txt+="""<plot title="PCG convergence%s" cmd="data=%s.pcg.tolerance" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="PCG alpha%s" cmd="data=%s.pcg.alphaHist" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="PCG beta%s" cmd="data=%s.pcg.betaHist" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        if self.parent.has_key("atmos"):
            #select this and run the simulation.  After an appropriate number of iterations, deselect.
            #This selects phaseCovariance, and sets atmos cal_source to 0.
            txt+="""<plot title="Start monte-phase covariance%s" when="cmd" ret="feedback" texttype="1" wintype="mainwindow">\n<cmd>\nfeedback=1-%s.control['phaseCovariance']\n%s.control['phaseCovariance']=feedback\nif feedback==0:\n %s.finishMontePhaseCovariance()\nfor obj in infAtmosList:\n obj.control['cal_source']=0\n obj.control['fullPupil']=feedback\n obj.control['removePiston']=1-feedback\n</cmd>\nbutton=feedback\nif feedback==0:\n msg='Phase covariance matrix computed'\nelse:\n msg='Computing phase covariance...'\nfeedback=msg</plot>\n"""%(id,objname,objname,objname)
            #Select this and run for an appropriate number of iterations, then deselect.
            txt+="""<plot title="Compute monte-noise covariance%s" when="cmd" ret="feedback" texttype="1" wintype="mainwindow">\n<cmd>\nfeedback=1-%s.control['covariance']\n%s.control['covariance']=feedback\nfor obj in infAtmosList:\n obj.control['cal_source']=feedback\nfor obj in dmList:\n obj.control['dm_update']=1\n%s.control['zero_dm']=feedback\n%s.control['close_dm']=1-feedback\n</cmd>\nif feedback==0:\n msg='Noise covariance matrix computed'\n button=0\nelse:\n msg='Computing covariance matrix'\n button=1\nfeedback=msg\n</plot>\n"""%(id,objname,objname,objname,objname)
            #and then once computed, we need to actually use them to compute the reconstructor...
            txt+="""<plot title="Compute monteCarlo reconstructor%s" when="cmd" cmd="%s.makeMonteRecon()"/>"""%(id,objname)
            txt+="""<plot title="Plot MC phase covariance%s" cmd="data=%s.montePhaseCovMat" ret="data" type="pylab" when="cmd"/>\n"""%(id,objname)
            txt+="""<plot title="Plot MC noise covariance%s" cmd="data=%s.subapVariance['cent']" ret="data" type="pylab" when="cmd"/>\n"""%(id,objname)
            txt+="""<plot title="Save MC noise covariance%s" cmd="import util.FITS;util.FITS.Write(%s.subapVariance['cent'],'subapVariance.fits');data='Saved'" ret="data" type="texttype" wintype="mainwindow" when="cmd"/>\n"""%(id,objname)
            txt+="""<plot title="Load MC noise covariance%s" cmd="import util.FITS;%s.subapVariance['cent']=util.FITS.Read('subapVariance.fits')[1];data='Loaded'" ret="data" type="texttype" wintype="mainwindow" when="cmd"/>\n"""%(id,objname)
            txt+="""<plot title="Load MC phase covariance%s" cmd="import util.FITS;%s.montePhaseCovMat=util.FITS.Read('montePhaseCovMat.fits')[1];data='Loaded'" ret="data" type="texttype" wintype="mainwindow" when="cmd"/>\n"""%(id,objname)
            txt+="""<plot title="Save reconmx.fits%s" cmd="import util.FITS;util.FITS.Write(%s.reconmx,'reconmx.fits');data='Saved'" ret="data" type="texttype" wintype="mainwindow" when="cmd"/>\n"""%(id,objname)
            txt+="""<plot title="Load reconmx.fits%s" cmd="import util.FITS;%s.reconmx=util.FITS.Read('reconmx.fits')[1];data='Loaded'" ret="data" type="texttype" wintype="mainwindow" when="cmd"/>\n"""%(id,objname)
            
            
        
        txt+="""<plot title="%s Output data%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt"/>\n"""%(self.objID,id,objname)

        return txt
    
if __name__=="__main__":
    print "Not testing xinterp_recon"
