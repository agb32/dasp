"""A module to compute the centroids of an input phase, as if this phase has been put into a SHS.
This can be used by the aosim simulation, or as standalone code.  When using
as stand alone code, typically you would create a centroid object, and then 
call the "run" method, passing in your phase.  Currently, using as stand alone
code means that the FPGAs cannot be easily used - though you can probably
change this easily enough if required.
"""
import util.flip,cmod.binimg,cmod.imgnoise
import numpy,numpy.random,numpy.fft
import util.arrayFromArray
import util.poisson,time,os
import util.dist
import util.centcmod
import util.correlation
haveFPGA=1
try:
    import fpga
except:
    print "FPGA module not installed."
    haveFPGA=0
haveCell=1
try:
    import util.centcell
except:
    print "cell centcell/centppu module not installed"
    haveCell=0


def pxlToRadTiltTheoretical(nsubx,phasesize,nimg,wfslam,telDiam):
    """compute the centroid pixel to radians of tilt scaling factor theoretically.  telDiam is in m, wfslam in nm.
    The logic behind this is as follows:
    Consider a collimated beam incident on a lens with focal length f.  The offset of the focussed
    spot in the image plane, x, is then x=f tan(t) where t is the angle of the incident wavefront.
    So, the question then is what is f in the simulation?  For each subaperture, there are n_p
    phase pixels, each with width z_p.  Similarly there are n_d detector pixels with width z_d.
    We note that n_p * z_p = n_d * z_d.
    Now, the airy disk created by the simulation is 2.44 * oversamplefactor pixels in diameter,
    where oversamplefactor is n_d/n_p.  The diameter of an airy spot is 2.44 lambda f/d = 2.44*l*f/(n_p*z_p).
    However, this is equal to the number of detector pixels wide, times the width of a detector
    pixel, i.e. 2.44*(n_d/n_p)*z_d.
    (Each pixel is lambda/d/oversamplefactor wide).  
    Equating these two expressions and cancelling terms therefore gives f = n_p * z_p^2/l.
    Therefore:
    t(x)=arctan(l * x / (n_p * z_p^2))

    So, multiply the centroids (in pixels) by this value to get the phase slopes in radians.
    When doing this, use x=1 pixel, ie x=telDiam/nsubx/nimg
    """
    npup=phasesize*nsubx
    x=telDiam/nsubx/nimg
    return numpy.arctan2((wfslam*1e-9*x),(phasesize*(telDiam/npup)**2.))#this is the radians of tilt that would give a 1 pixel shift in centroid algorithm.

def pxlToRadPistonTheoretical(nsubx,phasesize,nimg,wfslam,telDiam):
    """compute the centroid pixel to radians of piston scaling factor theoretically.  Here, the radians of piston are the wavefront that gives a certain tilt across a subaperture.
    """
    tilt=pxlToRadTiltTheoretical(nsubx,phasesize,nimg,wfslam,telDiam)#get radians of tilt that give a centroid shift of 1 pixel.
    #Then convert this tilt into a piston...
    x=numpy.tan(tilt)*(telDiam/nsubx)#piston in m...
    return x*2*numpy.pi/(wfslam*1e-9)
    
    
def wfs_sig(bandwidth,thruput,rateFactor,telDiam,nsubx,integTime,magnitude):
    """wfs signal per exposure
    bandwidth in angstroms (optical bandwidth)
    thruput including DQE etc
    rateFactor Photons/cm^2/s/Angstrom
    telDiam in m
    """
    pupilArea=(telDiam*100./nsubx)**2.      # Pupil area/cm^2
    sig=rateFactor*thruput*bandwidth*integTime*pupilArea/(2.5**magnitude) # Detected image photons
    #print 'WFS  photons/subap/integration: ',sig
    return sig

class centroid:
    """This class can be used as a library, to do wfs/centroid calculations.
    It is also used by the AO simulation, wfscent.py module.  So, any
    changes made here, will affect the simulation.
    This module can use the FPGAs.
    When used in the simulation framework, it enables resource sharing
    (memory and FPGA shared by more than one wfscent algorithm calculator).
    """
    def __init__(self,nsubx,pup=None,oversamplefactor=1,readnoise=0.,readbg=0.,addPoisson=0,noiseFloor=0.,binfactor=1,sig=1.,skybrightness=0.,warnOverflow=None,atmosPhaseType="phaseonly",fpDataType=numpy.float32,useFPGA=0,waitFPGA=0,waitFPGATime=0.,phasesize=None,fftsize=None,clipsize=None,nimg=None,ncen=None,tstep=0.05,integtime=0.05,latency=0.,wfs_minarea=0.5,spotpsf=None,centroidPower=None,opticalBinning=0,useCell=0,waitCell=1,usecmod=1,subtractTipTilt=0,magicCentroiding=0,linearSteps=None,stepRangeFrac=1.,phaseMultiplier=1,centWeight=None,correlationCentroiding=0,corrThresh=0.,corrPattern=None,threshType=0,imageOnly=0,calNCoeff=0,useBrightest=0):
        """
        Variables are:
         - sig: is the number of photons per phase pixel if pupfn is specified, or is the number
                of photons per subap if not (ie same...)  Or is a 2D array, with a value for each
                subap.  If this is just a number, it will be scaled by the number of used phase
                pixels for each subap.  If an array, assumes that this scaling has already been
                done.  Use an array version for eg LGS spot elongation.
         - nsubx: number of subaps in 1 direction
         - pupfn: pupil function array (pupil mask) or util.tel.Pupil instance (if using PS3).
         - oversamplefactor: scaling to expand phase by when doing FFTs.  Ignored if fftsize!=None.
         - readnoise: ccd readnoise
         - readbg: mean added signal from CCD.
         - addPoisson: whether to add photon shot noise
         - noiseFloor: Floor to apply to images.
         - binfactor: Bin factor to apply to get from FFT'd data to image.  Ignored if nimg!=None.
         - sig: signal
         - skybrightness: photons to add per pixel due to sky background
         - warnOverflow: will warn if the CCD is about to overflow (warnOverflow holds the saturation value)
         - atmosPhaseType: Type of phase from atmos module, eg phaseonly
         - fpDataType: data type used for calculations.
         - useFPGA: Whether FPGA should be used - shouldn't change.
         - waitFPGA: whether to wait for FPGA to finish before returning from calc.
         - waitFPGATime: time to wait before polling begins.
         - phasesize: wfs_n - number of phase pixels per subap
         - fftsize: wfs_nfft - number of pixels to use for FFT (zero pad phasesize...)
         - nimg: Number of pixels wide the subap image should be.
         - ncen: Number of subap image pixels used to compute centroid (eg could leave a dead ring around edge etc).
         - tstep: time per iteration
         - integtime: total integration time
         - latency: latency between readout starting and centroids being returned.
         - wfs_minarea: minimum fraction of subap receiving phase to allow it to be counted.
         - spotpsf: array (2 or 4d) of spot pattern PSFs.
         - centroidPower: None, or value to raise image to before centroiding.
         - opticalBinning: whether cylindrical lenslets are used, to do binning in 1 direction...
         - magicCentroiding: whether to measure phase slopes directly, or use SHS
         - linearSteps: None, or the number of steps to be used to try to use SH centroiding, but with a linear response (useful for openloop).  If None, linearisation is not carried out.
         - stepRangeFrac: The fraction of a subap over which to calibrate using linearSteps... default 1 for all of a subap.
         - phaseMultiplier can be non-unity if phase is at a different wavelength than you want to wfs at.
         - centWeight - weighting factor (array, 2 or 4D) for weight CoG centroiding.
         - correlationCentroiding - whether correlation centroiding is used.
         - corrThresh - subtraction threshold if using correlation cnetroiding
         - corrPattern - e.g. spot PSFs.
         - threshType - the threshold type for removing CCD read noise:
           If ==0, where(ccd<thresh:0:ccd-thresh)
           If ==1, where(ccd<thresh:0,ccd)
         - imageOnly - usually 0, but can be 1 or 2 if want to compute the CCD image only.
         - calNCoeff - used if linearSteps!=0, the number of coeffencients to use in polynomial fit (if this is zero, an interpolation routine is used instead)
         - useBrightest - int or array, if want to use brightest pixels algorithm.
        """
        self.nsubx=nsubx
        self.warnOverflow=warnOverflow
        self.sig=sig
        self.timing=0
        self.skybrightness=skybrightness
        self.pup=pup
        if type(pup)==numpy.ndarray or type(pup)==type(None):
            self.pupfn=pup
        else:
            self.pupfn=pup.fn
        self.oversamplefactor=oversamplefactor
        self.readnoise=readnoise
        self.readbg=readbg
        self.threshType=threshType
        self.imageOnly=imageOnly
        self.useBrightest=useBrightest
        self.binfactor=binfactor
        self.atmosPhaseType=atmosPhaseType
        self.fpDataType=fpDataType
        self.useFPGA=useFPGA#shouldn't change
        self.waitFPGA=waitFPGA
        self.waitFPGATime=waitFPGATime
        self.useCell=useCell#shouldn't change
        self.waitCell=waitCell
        self.phasesize=phasesize
        self.fftsize=fftsize
        self.clipsize=clipsize
        self.nimg=nimg
        self.ncen=ncen
        if self.clipsize==None:
            self.clipsize=self.fftsize
        if self.fftsize!=None and self.phasesize!=None:
            self.oversamplefactor=float(self.fftsize)/self.phasesize
        if self.clipsize!=None and self.nimg!=None:
            self.binfactor=float(self.clipsize)/self.nimg
        self.tstep=tstep
        self.integtime=integtime
        self.latency=latency
        self.nIntegrations=int(numpy.ceil(self.integtime/self.tstep))
        self.wfs_minarea=wfs_minarea
        self.psf=spotpsf#eg createAiryDisc(self.fftsize,self.fftsize/2,0.5,0.5)
        self.centroidPower=centroidPower
        self.opticalBinning=opticalBinning
        self.centWeight=centWeight

        #stuff for correlation centroiding...
        self.correlationCentroiding=correlationCentroiding
        self.corrThresh=corrThresh
        self.corrPatternUser=corrPattern
    
        if correlationCentroiding:
            self.corrimg=numpy.zeros((nsubx,nsubx,nimg,nimg),numpy.float32)
            self.corrPattern=util.correlation.transformPSF(self.corrPatternUser)
        else:
            self.corrPattern=None
            self.corrimg=None
        self.refCents=None

        self.convFactor=1.#convert from radians to centroids.
        self.printmax=0
        self.texp=0.#current exposure time.
        self.lastCalSource=0#value of calSource in last iteration (do we need to reload FPGA values?)
        self.addPoisson=addPoisson
        self.noiseFloor=noiseFloor
        self.usecmod=usecmod#should we use the cmodule version?
        self.subtractTipTilt=subtractTipTilt
        self.phasesize_v=self.phasesize#if using cell, this gets increased to a vectorised version (see below).
        if type(self.pupfn)!=type(None) and oversamplefactor!=None:
            self.subarea=numpy.ones((nsubx,nsubx),numpy.float64)
            n=self.pupfn.shape[0]/nsubx
            for i in range(nsubx):
                for j in range(nsubx):
                    self.subarea[i,j]=numpy.sum(numpy.sum(self.pupfn[i*n:(i+1)*n,j*n:(j+1)*n]))
        else:
            self.subarea=None
        self.wfsn=8#this is changed when phase is input...
        if self.nsubx!=1 and oversamplefactor!=None and binfactor!=None:
            self.computePxlToRad(self.wfsn)
        self.canUseFPGA=0
        self.canUseCell=0
        if self.useFPGA:
            self.sigFPGA=int(self.sig*2**7/(self.phasesize**2))&0x7ffffff#20.7 format.
            if int(self.sig*2**7)>0x7ffffff:#dont divide by phasesize^2 here.
                print "wfscent: Warning - SIG too large, will overflow in FPGA"
            self.fpgaSkybrightness=int(self.skybrightness*2**10/self.phasesize**2)
            self.fpga_readbg=int((self.readbg-self.readnoise*127/26.11)*256)
            self.fpga_readnoise=int(self.readnoise/26.11*256)
            self.testFPGAUsage()
            self.usePupil=1
            self.symtype=0
            npxls=(self.nsubx*self.phasesize)**2
            if npxls>4*1024*1024*8*4:
                print "Warning: can't fit pupil function into FPGA.  Assuming all pixels needed."
                self.usePupil=0
            elif npxls>4*1024*1024*8*2:
                print "Using 2 fold symmetry for pupil function"
                self.symtype=2
            elif npxls>4*1024*1024*8:
                print "Using 1 fold symmetry for pupil function"
                self.symtype=1
        elif self.useCell:
            self.canUseCell=haveCell
            self.phasesize_v=(self.phasesize+3)&~3#vectorised version for SPUs.
        self.magicCentroiding=magicCentroiding
        self.linearSteps=linearSteps
        self.calNCoeff=calNCoeff
        self.stepRangeFrac=stepRangeFrac
        self.phaseMultiplier=phaseMultiplier
        if magicCentroiding:
            self.magicSlopes=None#self.magicSHSlopes()

        self.subflag=numpy.zeros((self.nsubx,self.nsubx),numpy.int8)
        if self.subflag.itemsize==8:
            print "WARNING: untested with 8 byte longs...(wfs)"
        self.subarea=numpy.zeros((self.nsubx,self.nsubx),numpy.float64)
        n=self.phasesize
        #if self.pupfn==None:
        #    pfn=numpy.ones((1,1),self.fpDataType)
        #else:
        pfn=self.pupfn.astype(self.fpDataType)
        indices=[]
        for i in xrange(self.nsubx):        
            for j in xrange(self.nsubx):
                self.subarea[i,j]=pfn[i*n:(i+1)*n,j*n:(j+1)*n].sum()    # Get pupil fn over subaps
                if(self.subarea[i,j]>(self.wfs_minarea*n*n)):# Flag vignetted subaps 
                    self.subflag[i,j]=1
                    indices.append((i*self.nsubx+j)*2)
                    indices.append((i*self.nsubx+j)*2+1)
        self.indices=numpy.array(indices,dtype=numpy.int32)
        self.nsubaps=self.subflag.sum()
        #print "Created centroid object"

    def easy(self,nthreads=2,calsource=1):
        """Prepare for single use..., eg from a python commandline.
        Assumes that oversamplefactor=None, binfactor=None and that phasesize etc specified.
        """
        self.initMem(0)
        self.finishInit()
        self.initialiseCmod(nthreads,calsource)
        self.takeReference({'cal_source':1})
        self.outputData[:]=0
        #Then, put your data into self.reorderedPhs... 
        #Now, can use self.runCalc({'cal_source':0/1})

    def initMem(self,useFPGA,fpgaarr=None,shareReorderedPhs=0,subimgMem=None,bimgMem=None,pupsubMem=None,reorderedPhsMem=None,outputDataMem=None,useCell=0):
        """initialise memory banks - useful if resource sharing is used in the simulation.
        Not needed if not using simulation framework.
        if useFPGA is set, tells us that we're using the FPGA... (only used
        if allocating the memories here...)
        if useCell is set, tells us that we're using the cell...
        if all resource sharing objects expect to do the wfs/cent calc every iteration, then they can share reorderedPhs.  Otherwise, they must each have their own version.
        """
        #print "centroid - initMem"
        nsubx=self.nsubx
        phasesize=self.phasesize
        fftsize=self.fftsize
        nIntegrations=self.nIntegrations
        self.fittedSubaps=None
        self.shareReorderedPhs=shareReorderedPhs
        phasesize_v=self.phasesize_v
        if useFPGA:
            if self.atmosPhaseType=="phaseonly":
                atmosfactor=1#phaseonly...
            else:
                raise Exception("centroid: atmosphasetype...")
            if type(fpgaarr)==type(None):
                raise Exception("centroid: initMem called for FPGA use without passing the fpga accessible array.")
            self.fpgaarr=fpgaarr
            fpgaarrsize=fpgaarr.shape[0]
            self.fpgaarrsize=fpgaarrsize
            memNeeded=self.nsubx*self.nsubx*self.nIntegrations*self.phasesize*self.phasesize*4*atmosfactor+self.nsubx*self.nsubx*2*4
            nsubx=self.nsubx
            nIntegrations=self.nIntegrations
            phasesize=self.phasesize
            if memNeeded>1024*1024*1024:
                #cant all fit in fpga memory bank at once.
                self.doPartialFPGA=1
                #now see how many subaps can fit at once...
                self.fittedSubaps=1024*1024*1024/(nIntegrations*phasesize*phasesize*4*atmosfactor+2*4)
                self.partialFull=int(nsubx*nsubx/self.fittedSubaps)
                #and the number of subaps left to do.
                self.partialLeftOver=nsubx*nsubx-self.partialFull*self.fittedSubaps
                self.waitLeftOver=self.partialLeftOver*nIntegrations*phasesize*phasesize*5e-9
                #temporary input and output arrays for FPGA to access.
                if atmosfactor==1:
                    self.fpgaInArr=util.arrayFromArray.arrayFromArray(fpgaarr,(self.fittedSubaps,nIntegrations,phasesize,phasesize),numpy.float32)
                #create the input and output to copy from and to...
                #if self.atmosPhaseType!="phaseonly":
                #    pass
                self.fpgaOutArr=util.arrayFromArray.arrayFromArray(fpgaarr,(self.fittedSubaps,2),numpy.float32)
                if shareReorderedPhs==0 or type(reorderedPhsMem)==type(None):
                    self.reorderedPhs=numpy.zeros((nsubx,nsubx,nIntegrations,phasesize,phasesize),numpy.float32)
                else:
                    self.reorderedPhs=util.arrayFromArray.arrayFromArray(reorderedPhsMem,(nsubx,nsubx,nIntegrations,phasesize,phasesize),numpy.float32)
                if type(outputDataMem)==type(None):
                    if self.imageOnly==0:
                        self.outputData=numpy.zeros((nsubx,nsubx,2),numpy.float32)       # Centroid arrays
                    elif self.imageOnly==1:
                        self.outputData=numpy.zeros((nsubx,nsubx,self.nimg,self.nimg),numpy.float32)
                    else:
                        self.outputData=numpy.zeros((nsubx*self.nimg,nsubx*self.nimg),numpy.float32)
                                   
                else:
                    if self.imageOnly==0:
                        self.outputData=util.arrayFromArray.arrayFromArray(outputDataMem,(nsubx,nsubx,2),numpy.float32)
                    elif self.imageOnly==1:
                        self.outputData=util.arrayFromArray.arrayFromArray(outputDataMem,(nsubx,nsubx,self.nimg,self.nimg),numpy.float32)
                    else:
                        self.outputData=util.arrayFromArray.arrayFromArray(outputDataMem,(nsubx*self.nimg,nsubx*self.nimg),numpy.float32)
                        
            else: #fpga can access whole array.
                self.doPartialFPGA=0
                #fpgaarr=fpga.mallocHostMem(fpid,(fpgaarrsize,),numpy.Int8)
                #if self.atmosPhaseType!="phaseonly":
                #    pass
                if shareReorderedPhs==0:
                    self.reorderedPhs=numpy.zeros((nsubx,nsubx,nIntegrations,phasesize,phasesize),numpy.float32)
                    self.fpgaInput=util.arrayFromArray.arrayFromArray(fpgaarr,(nsubx,nsubx,nIntegrations,phasesize,phasesize),numpy.float32)#copy phase to here, before running the fpga.
                else:#it has to come from the fpgaarr... 
                    self.reorderedPhs=util.arrayFromArray.arrayFromArray(fpgaarr,(nsubx,nsubx,nIntegrations,phasesize,phasesize),numpy.float32)#fpga can access this directly.
                self.outputData=util.arrayFromArray.arrayFromArray(fpgaarr[nsubx*nsubx*nIntegrations*phasesize*phasesize*4:,],(nsubx,nsubx,2),numpy.float32)
        else:#not using FPGA, so set up memory without it.
            if shareReorderedPhs==0 or type(reorderedPhsMem)==type(None):
                self.reorderedPhs=numpy.zeros((nsubx,nsubx,nIntegrations,phasesize,phasesize_v),numpy.float32)
            else:
                self.reorderedPhs=util.arrayFromArray.arrayFromArray(reorderedPhsMem,(nsubx,nsubx,nIntegrations,phasesize,phasesize_v),numpy.float32)
            if type(outputDataMem)==type(None):
                if self.imageOnly==0:
                    self.outputData=numpy.zeros((nsubx,nsubx,2),numpy.float32)       # Centroid arrays
                elif self.imageOnly==1:
                    self.outputData=numpy.zeros((nsubx,nsubx,self.nimg,self.nimg),numpy.float32)
                else:
                    self.outputData=numpy.zeros((nsubx*self.nimg,nsubx*self.nimg),numpy.float32)
                                   
            else:
                if self.imageOnly==0:
                    self.outputData=util.arrayFromArray.arrayFromArray(outputDataMem,(nsubx,nsubx,2),numpy.float32)
                elif self.imageOnly==1:
                    self.outputData=util.arrayFromArray.arrayFromArray(outputDataMem,(nsubx,nsubx,self.nimg,self.nimg),numpy.float32)
                else:
                    self.outputData=util.arrayFromArray.arrayFromArray(outputDataMem,(nsubx*self.nimg,nsubx*self.nimg),numpy.float32)
        if self.imageOnly==0:
            self.centx=self.outputData[:,:,0]
            self.centy=self.outputData[:,:,1]
        else:
            self.centx=None
            self.centy=None
        #self.outputData.savespace(1)
        if type(subimgMem)==type(None):
            # SH sub-images (high LL):
            self.subimg=numpy.zeros((self.nsubx,self.nsubx,self.fftsize,self.fftsize),numpy.float64)
        else:
            self.subimg=util.arrayFromArray.arrayFromArray(subimgMem,(self.nsubx,self.nsubx,self.fftsize,self.fftsize),numpy.float64)
        if type(bimgMem)==type(None):
            self.bimg=numpy.zeros((self.nsubx,self.nsubx,self.nimg,self.nimg),numpy.float64)
        else:
            self.bimg=util.arrayFromArray.arrayFromArray(bimgMem,(self.nsubx,self.nsubx,self.nimg,self.nimg),numpy.float64)
        #self.bimg.savespace(1)
        if self.usecmod:
            self.cmodbimg=util.arrayFromArray.arrayFromArray(self.bimg,(self.nsubx,self.nsubx,self.nimg,self.nimg),numpy.float32)
        else:
            self.cmodbimg=None
##         if type(shimgMem)==type(None):
##             self.shimg=numpy.zeros((self.nsubx*self.nimg,self.nsubx*self.nimg),numpy.float64)# Tessalated SH image for display
##         else:
##             self.shimg=util.arrayFromArray.arrayFromArray(shimgMem,(self.nsubx*self.nimg,self.nsubx*self.nimg),numpy.float64)
        if type(pupsubMem)==type(None):
            self.pupsub=numpy.zeros((self.nsubx,self.nsubx,self.phasesize,self.phasesize),self.fpDataType)
        else:
            self.pupsub=util.arrayFromArray.arrayFromArray(pupsubMem,(self.nsubx,self.nsubx,self.phasesize,self.phasesize),self.fpDataType)
        #print "centroid - initMem done"
    def initialiseCell(self,nspu=6,calsource=0,showCCDImg=0,allCents=1,cellseed=1):
        if self.canUseCell:
            self.cellObj=util.centcell.centcell(self.fftsize,self.nsubx,self.nimg,self.phasesize,
                                                self.ncen,self.nIntegrations,self.reorderdPhs,self.psf,self.pup,
                                                self.outputData,calsource,self.sig,self.readnoise,
                                                readbg=self.readbg,noisefloor=self.noiseFloor,
                                                seed=cellseed,minarea=self.wfs_minarea,
                                                skybrightness=self.skybrightness,allCents=allCents,
                                                shimg=None,nspu=nspu)
            self.cellObj.showCCDImg=showCCDImg
            self.cellObj.initialise()

    def initialiseCmod(self,nthreads=8,calsource=0,seed=1):
        self.nthreads=nthreads
        if self.usecmod:
            print "initialising cmod, nthreads = {0}".format(nthreads)
            sig=self.sig
            if type(self.sig)==numpy.ndarray:
                sig=self.sig.ravel()
            #temporary til I get it done properly...
            self.centcmod=util.centcmod.centcmod(nthreads,self.nsubx,self.ncen,self.fftsize,self.clipsize,
                                                 self.nimg,self.phasesize,self.readnoise,self.readbg,
                                                 self.addPoisson,self.noiseFloor,sig,self.skybrightness,
                                                 calsource,self.centroidPower,self.nIntegrations,seed,
                                                 self.reorderedPhs,self.pup,self.psf,self.outputData,
                                                 self.cmodbimg,self.wfs_minarea,self.opticalBinning,
                                                 self.centWeight,self.correlationCentroiding,
                                                 self.corrThresh,self.corrPattern,self.corrimg,
                                                 self.threshType,self.imageOnly,self.useBrightest)
            #print "initialised cmod - done"
        else:
            self.centcmod=None
            
    def initialiseFPGA(self,fpid=None,ignoreFailure=0,fpgaInfo=None,fpgaBitFile=None):
        """Load the bin file.  Not needed if not using FPGAs."""
        self.fpid=fpid
        self.fpgaInfo=fpgaInfo#store info about what is currently loaded in fpga - eg fft size, pupil map etc.  This should be shared by all centroid objects using the FPGA, and is an instance of fpgaCentStateInformation class.
        self.fpgaBinaryLoaded=0
        if self.canUseFPGA:
            if self.fpid==None:
                try:
                    self.fpid=fpga.open(reset=0,start=0)
                    self.fpgaBinaryLoaded=1
                except:
                    if ignoreFailure:
                        self.fpgaBinaryLoaded=0
                        self.fpid=None
                        print "Warning: wfscent - failed to initialise FPGA"
                    else:
                        raise
                if self.fpgaBinaryLoaded:
                    fpga.load(self.fpid,fpgaBitFile)
                    if os.environ["HOSTNAME"]=="n1-c437":
                        print "WARNING: wfscent FPGA module may not work correctly on node1... possibly something not quite right with the FPGA."
                    fpga.reset(self.fpid)
                    time.sleep(0.001)
                    fpga.start(self.fpid)
                    time.sleep(0.001)
                    fpga.writeReg(self.fpid,0x2,6)#stop the fpga pipeline
                    self.fpgaInfo=fpgaCentStateInformation()
            else:
                self.fpgaBinaryLoaded=1
        if self.fpgaBinaryLoaded==0:
            self.fpid=None
            self.fpgaInfo=None
            self.canUseFPGA=0
        return self.fpid,self.fpgaInfo

    
    def setupFPGAArray(self,fpid=None,fpgaarrsize=None):
        """allocate FPGA memory buffer.  This assumes that phasetype is "phaseonly".
        """
        #fpgaarrsize,doPartialFPGA,fittedSubaps,partialLeftOver,partialFull,waitLeftOver,fpgaarr,fpgaInArr,reorderedPhs,outputData
        #if self.atmosPhaseType=="phaseonly":
        atmosfactor=1#phaseonly...
        #else:
        #    raise Exception("centroid: atmosphasetype...")
        nsubx=self.nsubx
        nIntegrations=self.nIntegrations
        phasesize=self.phasesize
        if fpgaarrsize==None:
            fpgaarrsize=nsubx*nsubx*nIntegrations*phasesize*phasesize*4*atmosfactor+nsubx*nsubx*2*4
        if fpid==None:
            fpid=self.fpid
        if fpgaarrsize<4*1024*1024:
            fpgaarrsize=4*1024*1024#needs 4MB array for loading QDR memory.
        if fpgaarrsize>1024*1024*1024:
            print "Warning: Pupil pixel size is too large for single FPGA array - will use multiple arrays, but speed will be reduced (FPGA can access only a 1GB buffer)."
            fpgaarrsize=1024*1024*1024#might not all be used, but then we don't know that!  However, it is likely that a whole number of subaps will fit exactly since they are usually powers of two.
        
        fpgaarr=fpga.mallocHostMem(fpid,(fpgaarrsize,),numpy.int8)
        return fpgaarr
        




    def finishInit(self):
        """other initialisations to be carried out after initMem has been called.
        Not needed if not using simulation framework"""
        #print "centroid - finishInit"
        self.xtiltfn=((numpy.fromfunction(self.tilt,(self.phasesize,self.phasesize))-float(self.phasesize)/2.+0.5)/float(self.phasesize)).astype(self.fpDataType)# subap tilt fn
        self.ytiltfn=numpy.transpose(self.xtiltfn)
        #if self.fpDataType==numpy.float64:
        #    self.fftwPlan=cmod.mkimg.setup(self.subimg[0][0])                   # Setup imaging FFTs for subaps
        #else:
        #    self.fftwPlan=cmod.mkimgfloat.setup(self.subimg[0][0])
        self.tmpImg=numpy.zeros((self.fftsize,self.fftsize),self.fpDataType)
        n=self.phasesize
        pfn=self.pupfn.astype(self.fpDataType)
        for i in xrange(self.nsubx):        
            for j in xrange(self.nsubx):
                self.pupsub[i][j]=pfn[i*n:(i+1)*n,j*n:(j+1)*n]    # Get pupil fn over subaps

        self.cenmask=numpy.zeros((self.nimg,self.nimg),numpy.float32)             # Define centroiding mask
        self.cenmask[self.nimg/2-self.ncen/2:self.nimg/2+self.ncen/2,self.nimg/2-self.ncen/2:self.nimg/2+self.ncen/2]=1.
        self.tilt_indx = (numpy.array(range(self.nimg),numpy.float64))-float(self.nimg/2)+0.5#Index fns for centroiding
        #print "centroid - finishInit done"

    def closeCell(self):
        if self.canUseCell:
            self.cellObj.close()

    def runCalc(self,control):
        doref=1
        if self.phaseMultiplier!=1:
            self.reorderedPhs*=self.phaseMultiplier
        if control.get("useFPGA",0) and self.canUseFPGA:
            # use the FPGA - note that you might get a non-zero centroid value for parts of
            # the array which are masked off simply because of the ccd readout noise.  The 
            # software version handles this incorrectly.
            # Check whether registers are still valid for this object, and if not, change them so that they are:
            self.setFPGARegs(control["cal_source"])
            self.runFPGA()
        elif control.get("useCell",0) and self.canUseCell:
            self.runCell(control["cal_source"])
        elif control.get("useCmod",1):
            self.runCmod(control["cal_source"])
            # no calibration done, or done in c, so ref can be done by c:
            if self.linearSteps==None or self.psf!=None or self.correlationCentroiding!=None or self.calNCoeff!=0:
                doref=0#ref subtraction is done in c code...
        else:
            # use software version
            # t=time.time()
            # Create the images
            self.runPy(control["cal_source"])
        if self.linearSteps!=None:
            self.applyCalibration()
        if doref:
            if self.refCents!=None:
                self.outputData-=self.refCents

    def runCell(self,calsource):
        """Tell the cell to perform computations."""
        if self.magicCentroiding:
            self.magicShackHartmann()
            return
        if self.canUseCell:
            self.cellObj.setCalSource(calsource)
            # If waitCell==0, will need to call self.cellObj.waitForCents() at some later time:
            self.cellObj.startProcessing(block=self.waitCell)

    def runPy(self,calsource):
        """run the python version"""
        if self.magicCentroiding:
            self.magicShackHartmann()
            return
        self.createSHImgs()
        self.tidyImage(calsource)
        self.calc_cents(calsource)

    def runCmod(self,calsource):
        """run the c version"""
        if self.magicCentroiding:
            self.magicShackHartmann()
            return
        self.centcmod.run(calsource)
        if self.imageOnly==0:
            pass
        elif self.imageOnly==1:
            #copy image into outputdata
            self.outputData[:]=self.cmodbimg
        else:
            #copy image as image into outputdata
            for i in xrange(self.nsubx):
                for j in xrange(self.nsubx):
                    self.outputData[i*self.nimg:(i+1)*self.nimg,j*self.nimg:(j+1)*self.nimg]=self.cmodbimg[i,j]

    def closeCmod(self):
        self.centcmod.free()
        self.centcmod=None

    def runFPGA(self):
        """Tell the FPGA where the data is..."""
        if self.magicCentroiding:
            self.magicShackHartmann()
            return
        fpid=self.fpid
        #now, we set the FPGA going (after copying data if necessary).
        t0=time.time()
        #print "runfpga"
        if self.doPartialFPGA:
            #array is too large to DMA all at once to FPGA, so do in parts.
            #calculate number of times a full array is needed...
            if self.atmosPhaseType=="phaseonly":
                reordered=util.arrayFromArray.arrayFromArray(self.reorderedPhs,
                          (self.nsubx*nsubx,self.nIntegrations,self.phasesize,self.phasesize),
                                                             numpy.float32)
            else:
                raise Exception("not phaseonly")
            output=util.arrayFromArray.arrayFromArray(self.outputData,(self.nsubx*self.nsubx,2),numpy.float32)
            fpga.writeAddr(fpid,self.fpgaInArr,1)#input address
            fpga.writeReg(fpid,self.fittedSubaps*self.nIntegrations*self.phasesize*self.phasesize*4/8,2)#size (quad words)
            fpga.writeAddr(fpid,self.fpgaOutArr,3)#output address
            fpga.writeReg(fpid,self.fittedSubaps,4)#size to write (in bytes).
            for i in range(self.partialFull):
                #copy memory into FPGA buffer
                self.fpgaInArr[:,]=reordered[i*self.fittedSubaps:(i+1)*self.fittedSubaps]
                fpga.writeReg(fpid,0x2,6)#reinitialise
                fpga.writeReg(fpid,1,6)#set it going.
                time.sleep(self.waitFPGATime/self.partialFull)#wait for it to complete (or almost)
                while fpga.readReg(fpid,5)!=7:#wait for reading to complete by checking register.
                    pass
                #copy centroids to the output...
                output[i*self.fittedSubaps:(i+1)*self.fittedSubaps]=self.fpgaOutArr
            if self.partialLeftOver>0:#copy the last bit...
                self.fpgaInArr[:self.partialLeftOver]=reordered[self.partialFull*self.fittedSubaps:self.partialFull*self.fittedSubaps+self.partialLeftOver]
                fpga.writeReg(fpid,0x2,6)#reinitialise
                fpga.writeReg(fpid,self.partialLeftOver*self.nIntegrations*self.phasesize*self.phasesize*4/8,2)#read siz
                fpga.writeReg(fpid,self.partialLeftOver,4)#size to write
                fpga.writeReg(fpid,1,6)#set it going
                time.sleep(self.waitLeftOver)
                while fpga.readReg(fpid,5)!=7:#wait til finished
                    pass
                #and copy centroids to the output array.
                output[self.partialFull*self.fittedSubaps:self.partialFull*self.fittedSubaps+self.partialLeftOver]=self.fpgaOutArr[:self.partialLeftOver]

                #now reset the registers for next time...
                fpga.writeReg(fpid,self.fittedSubaps*self.nIntegrations*self.phasesize*self.phasesize*4/8,2)#read siz
                fpga.writeReg(fpid,self.fittedSubaps,4)#size to write

                      
        else:
            #all subaps at once...
            if self.shareReorderedPhs:#this must be the only object using it...
                pass#already in the fpga array
            else:#copy to fpga array.
                self.fpgaInput[:,]=self.reorderedPhs
            fpga.writeReg(fpid,0x2,6)#reinitialise
            fpga.writeAddr(fpid,self.fpgaarr,1)#input address
            fpga.writeReg(fpid,self.nsubx*self.nsubx*self.nIntegrations*self.phasesize*self.phasesize*4/8,2)#size (quad words)
            fpga.writeAddr(fpid,self.outputData,3)#output address
            fpga.writeReg(fpid,self.nsubx*self.nsubx,4)#size to write (in bytes).
            #print "reading fpga reg %s"%hex(fpga.readReg(fpid,5))
            t0=time.time()
            fpga.writeReg(fpid,1,6)#set it going.
            if self.waitFPGA:
                if self.waitFPGATime>0:
                    time.sleep(self.waitFPGATime)#wait for it to complete (or almost).
                v=fpga.readReg(fpid,5)
                #print hex(v)
                while v!=7:#wait for reading to complete by checking register...
                    v=fpga.readReg(fpid,5)
                    #print hex(v)
                    pass
        if self.timing:
            t1=time.time()
            print "WFSCent time taken: %s"%str(t1-t0)
        #print "runfpgadone"
    def reorder(self,phs,pos):
        """Do a reodering of the phase buffer, so that it is in the form ready
        for the fpga.  Phases corresponding to subaps should be placed next to
        those for the same subap at a later time (greater pos).
        Also use with createSHImg() but not integrate().  If the FPGA has been initialised, this will place it in FPGA readable memory..."""
        #nsubx=    self.nsubx
        n=    self.phasesize
        typecode=self.reorderedPhs.dtype
        if self.atmosPhaseType=="phaseonly":
            for i in range(self.nsubx):
                for j in range(self.nsubx):
                    # Subap phase array
                    self.reorderedPhs[i,j,pos,:,:n]=phs[i*n:(i+1)*n,j*n:(j+1)*n].astype(typecode)
        elif self.atmosPhaseType=="phaseamp":
            for i in range(self.nsubx):
                for j in range(self.nsubx):
                    # Subap phase array
                    self.reorderedPhs[i,j,pos,:,:n,0]=phs[0,i*n:(i+1)*n,j*n:(j+1)*n].astype(typecode)#phase
                    self.reorderedPhs[i,j,pos,:,:n,1]=phs[1,i*n:(i+1)*n,j*n:(j+1)*n].astype(typecode)#amplitude
        else:
            raise Exception("atmosphasetype")

    def makeImage(self,phs,img,pup):
        tmpphs=numpy.zeros(img.shape,numpy.complex64)
        tmpphs[:phs.shape[0],:phs.shape[1]]=(pup*(numpy.cos(phs)+1j*numpy.sin(phs))).astype(numpy.complex64)
        tmp=numpy.fft.fft2(tmpphs)
        img[:,]=(tmp.real*tmp.real+tmp.imag*tmp.imag).astype(numpy.float32)

    def createSHImgs(self):
        """Do the FFTs etc and add integrations to create the powerspectra (high light level images).
        All this will eventually (we hope) be done in the FPGAs"""
        tmp=0.5*float(self.phasesize)/float(self.fftsize)*2.*numpy.pi
        self.subimg*=0.0                                                  # Zero the CCD
        #nsubx=self.nsubx
        for i in xrange(self.nsubx):
            for j in xrange(self.nsubx):
                if self.subflag[i][j]==1:
                    for k in xrange(self.nIntegrations):
                        if self.atmosPhaseType=="phaseonly":
                            phssub=(self.reorderedPhs[i,j,k,:,:self.phasesize]-tmp*self.xtiltfn-tmp*self.ytiltfn).astype(self.fpDataType)
                        elif self.atmosPhaseType=="phaseamp":
                            phssub=(self.reorderedPhs[i,j,k,:,:,0]-tmp*self.xtiltfn-tmp*self.ytiltfn).astype(self.fpDataType)
                        #now do the FFT: plan, phsin, imgout, pypupil
                        if self.fpDataType==numpy.float64:
                            if phssub.dtype!=self.fpDataType or self.tmpImg.dtype!=self.fpDataType or self.pupsub.dtype!=self.fpDataType:
                                print "ERROR with typecodes in wfs.py"
                                raise Exception("Error with typecodes in wfs.py mkimg")
                            if self.atmosPhaseType=="phaseonly":
                                #cmod.mkimg.mkimg(self.fftwPlan,phssub,self.tmpImg,self.pupsub[i][j])
                                self.makeImage(phssub,self.tmpImg,self.pupsub[i,j])
                            elif self.atmosPhaseType=="phaseamp":
                                #cmod.mkimg.mkimg(self.fftwPlan,phssub,self.tmpImg,(self.pupsub[i,j]*self.reorderedPhs[i,j,k,:,:,1]).astype(self.fpDataType))
                                self.makeImage(pupsub,self.tmpImg,(self.pupsub[i,j]*self.reorderedPhs[i,j,k,:,:,1]).astype(self.fpDataType))
                            elif self.atmosPhaseType=="realimag":
                                raise Exception("realimag")
                        else:
                            if phssub.dtype!=self.fpDataType or self.tmpImg.dtype!=self.fpDataType or self.pupsub.dtype!=self.fpDataType:
                                print "ERROR with typecodes in wfs.py",phssub.dtype,self.tmpImg.dtype,self.pupsub.dtype
                                raise Exception("Error with typecodes in wfs.py mkimg")
                            if self.atmosPhaseType=="phaseonly":
                                #cmod.mkimgfloat.mkimg(self.fftwPlan,phssub,self.tmpImg,self.pupsub[i][j])
                                self.makeImage(phssub,self.tmpImg,self.pupsub[i,j])
                            elif self.atmosPhaseType=="phaseamp":
                                #cmod.mkimgfloat.mkimg(self.fftwPlan,phssub,self.tmpImg,(self.pupsub[i,j]*self.reorderedPhs[i,j,k,:,:,1]).astype(self.fpDataType))
                                self.makeImage(phssub,self.tmpImg,(self.pupsub[i,j]*self.reorderedPhs[i,j,k,:,:,1]).astype(self.fpDataType))
                            elif self.atmosPhaseType=="realimag":
                                raise Exception("realimag")
                        self.subimg[i][j]+=self.tmpImg                    # Long exposure image

    def tidyImage(self,calSource):
        """Flip image (to correct for FFT), then bin the image up if
        oversampled, then add read and shot noise."""
        
        #self.shimg*=0.                                                # Reset...=zeros((nsubx*nimg,nsubx*nimg),float)
        readnoise=  self.readnoise
        mean=self.readbg
##         m=0
##         mn=1000000
        #nsubx=self.nsubx
        for i in range(self.nsubx):                                        # Loop over subaps
            for j in range(self.nsubx):
                if(self.subflag[i][j]==1):
                    bimg=self.bimg[i][j]
                    img=util.flip.fliparray(self.subimg[i][j])                          # Get subap image, flip it
                    #agb - test convolution... is this the right place?
                    if type(self.psf)!=type(None):
                        if len(self.psf.shape)==4:
                            ca=self.psf[i,j]
                        else:
                            ca=self.psf
                        img[:,]=self.conv(img,ca)

                    cmod.binimg.binimg(img,bimg)                              # Bin it up
                    totsig = numpy.sum(numpy.sum(bimg))
                    nphs=float(self.subarea[i,j])/self.phasesize**2#fraction of active phase pixels
                    if(totsig>0.):
                        if type(self.sig)==type(0.0):
                            bimg*=self.sig*nphs/totsig  #rescale the image.
                            # Note, used to be an error here - should also scale by the number of pixels that receive photons.  Fixed...
                        else:
                            bimg*=self.sig[i,j]/totsig
                        if self.atmosPhaseType!="phaseonly":
                            print "Scaling of SH image? Is it needed when have phase an amplitude for atmosphere"
                    #now add sky brightness scaled by the number of phase pixels that are in the aperture:
                    if not calSource:#self.control["cal_source"]:
                        bimg+=self.skybrightness*nphs
                        if self.opticalBinning:
                            s1=numpy.sum(bimg)/2#bin in 1 dimension and take half the light (beam splitter).
                            s2=numpy.sum(numpy.transpose(bimg))/2#bin in other dimension
                            bimg[0]=s1#and store back in bimg.
                            bimg[1]=s2
                            bimg[2:,]=0.#this isn't necessary but makes nicer image...
                            bimg=bimg[:2]
                        if totsig>0. or self.skybrightness*nphs>0.:
                            # Inject shot noise
                            cmod.imgnoise.shot(bimg,bimg)


                        # Generate random read noise :
                        bimg+=(numpy.random.normal(mean,readnoise,bimg.shape)+0.5).astype("i")#round to integer
                    #self.shimg[i*self.wfs_nimg:(i+1)*self.wfs_nimg,j*self.wfs_nimg:(j+1)*self.wfs_nimg]=bimg    # Tessalate up for WFS display

    def calc_cents(self,calSource):
        """Centroid calculation:
        Subtracts noise background
        Computes centroids
        No return value."""
        nfft=   self.fftsize
        nimg=self.nimg
        nsubx=  self.nsubx
        floor=  self.noiseFloor
        #read=  self.wfs_read
        indx=   self.tilt_indx
        #bimg=self.bimg
        self.outputData[:,]=0.                                      # Reset...
        #self.centx=zeros((nsubx,nsubx),float)
        #self.centy=zeros((nsubx,nsubx),float)
        #self.shimg*=0.#reset...=zeros((nsubx*nimg,nsubx*nimg),float)
        cenmask=self.cenmask
        if self.opticalBinning:
            cenmask=1
        for i in range(nsubx):                                  # Loop over subaps
            for j in range(nsubx):
                if(self.subflag[i][j]==1):
                    bimg=self.bimg[i,j]
                    if self.opticalBinning:
                        bimg=bimg[:2]
                    #threshold the SH images.
                    if not calSource:#self.control["cal_source"]:
                        cimg=numpy.where(bimg<self.noiseFloor,0,bimg-self.noiseFloor)*cenmask
                    else:
                        cimg=bimg*cenmask
                    if self.opticalBinning:#using a cylindrical lenslet array...
                        s1=numpy.sum(cimg[0])
                        if s1>0:
                            self.centx[i,j]=numpy.sum(cimg[0]*indx)/s1
                        else:
                            self.centx[i,j]=0.
                        s1=numpy.sum(cimg[1])
                        if s1>0:
                            self.centy[i,j]=numpy.sum(cimg[1]*indx)/s1
                        else:
                            self.centy[i,j]=0.
                    else:
                        totsig = numpy.sum(numpy.sum(cimg))
                        if(totsig==0.):
                            totsig=1.#division by zero...
                        # Centroid calculation
                        self.centx[i,j]=numpy.sum(numpy.sum(cimg,0)*indx)/totsig  
                        self.centy[i,j]=numpy.sum(numpy.sum(cimg,1)*indx)/totsig

    def testFPGAUsage(self):
        """Checks variables are suitable for FPGA use"""
        self.canUseFPGA=haveFPGA
        if self.atmosPhaseType!="phaseonly":
            print "WARNING: Cannot use FPGA - atmosPhaseType must be phaseonly"
            self.canUseFPGA=0
        if self.fpDataType!=numpy.float32:
            print "WARNING: Cannot use FPGA - fpDataType must be float32"
            self.canUseFPGA=0
        if self.fftsize not in [8,16,32]:
            print "WARNING: FPGA cannot use this FFT array size"
            self.canUseFPGA=0
        if self.nIntegrations<1 or self.nIntegrations>63:
            print "WARNING: Illegal number of integrations for FPGA - must be less than 64"
            self.canUseFPGA=0
        if self.nimg<2 or self.nimg>32:
            print "WARNING: Illegal pixel size for centroiding in FPGA (nimg) - must be 2-32"
            self.canUseFPGA=0
        if self.phasesize<1 or self.phasesize>32:
            print "WARNING: Illegal phase size for use in FPGA (phasesize) - must be <32"
            self.canUseFPGA=0
        if type(self.sig)!=type(0.0):
            print "WARNING: Signal cannot be array for use in FPGA."
            self.canUseFPGA=0
        if int(self.sig*2**7)>0x7ffffff:
            print "WARNING: Signal is too bright for use in FPGA - should be less than 0xfffff"
            self.canUseFPGA=0
        if self.skybrightness>0xffff:
            print "WARNING: Sky background is too bright for use in FPGA - should be less than 0xffff"
            self.canUseFPGA=0
        if self.noiseFloor>0xffff:
            print "WARNING: WFS floor (threshold) value is too high for use in FPGA - should be less than 0xffff"
            self.canUseFPGA=0
        if self.fpga_readbg>0xffffff or self.fpga_readnoise>0xffff:
            print "WARNING: CCD readout noise is too high for use in FPGA (mean or sigma)"
            self.canUseFPGA=0
        if self.nsubx<1 or self.nsubx>1023:
            print "WARNING: Number of subapertures is too large for use in FPGA (needs nsubx<1024)"
            self.caUseFPGA=0
        if type(self.psf)!=type(None):
            print "WARNING: FPGA cannot use psf (eg lgs spot elongation, or airy disc convolution)"
            self.canUseFPGA=0
        if self.opticalBinning:
            print "WARNING: FPGA cannot use optical binning"
            self.canUseFPGA=0
        return self.canUseFPGA

    def setFPGARegs(self,calsource):
        """Here we check that the registers in the FPGA are what they should be for this instance of centroid.  If not, we set the registers correctly.  If using resource sharing, this should be called everytime a new FPGA calc is required."""
        #print "setfpgaregs"
        if not self.canUseFPGA:
            return
        #set up memory...
        if self.fpgaInfo.calSource!=calsource:
            self.loadPoisson(calsource)
        if self.fpgaInfo.seedLoaded==0:
            self.loadSeed()
        if self.fpgaInfo.pupil is not self.pupfn:#different arrays
            self.loadPupil()
        #now initialise the registers...
        if self.fpgaInfo.nimg!=self.nimg or self.fpgaInfo.fftsize!=self.fftsize:
            self.loadBinCtrl()
        if type(self.fpgaInfo.cenmask)==type(None) or self.fpgaInfo.cenmask.shape!=self.cenmask.shape or numpy.sum(numpy.sum(self.fpgaInfo.cenmask==self.cenmask))!=self.nimg*self.nimg:
            self.loadCenMask()
        if self.fpgaInfo.doPartialFPGA!=self.doPartialFPGA or self.fittedSubaps!=self.fpgaInfo.fittedSubaps or self.nIntegrations!=self.fpgaInfo.nIntegrations or self.nsubx!=self.fpgaInfo.nsubx or self.phasesize!=self.fpgaInfo.phasesize:
            self.loadFPGARegAddr()

        if self.nIntegrations!=self.fpgaInfo.nIntegrations or self.nimg!=self.fpgaInfo.nimg or self.fftsize!=self.fpgaInfo.fftsize or self.phasesize!=self.fpgaInfo.phasesize or self.nsubx!=self.fpgaInfo.nsubx or self.symtype!=self.fpgaInfo.symtype:
            self.loadFPGADimData()
            
        if self.sigFPGA!=self.fpgaInfo.sigFPGA or calsource!=self.fpgaInfo.calSource or self.fpgaSkybrightness!=self.fpgaInfo.fpgaSkybrightness or self.noiseFloor!=self.fpgaInfo.noiseFloor or self.fpga_readbg!=self.fpgaInfo.fpga_readbg or self.fpga_readnoise!=self.fpgaInfo.fpga_readnoise:
            self.loadFPGASourceData(calsource)
        fpid=self.fpid
        fpga.writeReg(fpid,0x2,6)#stop pipe
        fpga.writeReg(fpid,64,6)#reset input/output fifos.
        fpga.writeReg(fpid,2,512)#set pipe to do WFSing.
        
        self.fpgaInfo.calSource=calsource
        self.fpgaInfo.seedLoaded=1#never needs reloading.
        self.fpgaInfo.pupil=self.pupfn
        self.fpgaInfo.nimg=self.nimg
        self.fpgaInfo.fftsize=self.fftsize
        self.fpgaInfo.cenmask=self.cenmask
        self.fpgaInfo.doPartialFPGA=self.doPartialFPGA
        self.fpgaInfo.fittedSubaps=self.fittedSubaps
        self.fpgaInfo.nIntegrations=self.nIntegrations
        self.fpgaInfo.nsubx=self.nsubx
        self.fpgaInfo.phasesize=self.phasesize
        self.fpgaInfo.sigFPGA=self.sigFPGA
        self.fpgaInfo.fpgaSkybrightness=self.fpgaSkybrightness
        self.fpgaInfo.noiseFloor=self.noiseFloor
        self.fpgaInfo.fpga_readbg=self.fpga_readbg
        self.fpgaInfo.fpga_readnoise=self.fpga_readnoise
        self.fpgaInfo.symtype=self.symtype

    def loadPupil(self):
        """Load the pupil map - note, using this overwrites the reorderedPhs array..."""
        print "Loading FPGA pupil map"
        #Now, place the pupil into a bit format that can be read by FPGA.
        usexsubap=self.nsubx
        useysubap=self.nsubx
        self.puparr=numpy.zeros((4*1024*1024,),numpy.uint8)
        if self.usePupil==0:#use every pixel.
            self.puparr[:,]=0xff
        else:
            if self.symtype==0:
                pupfn=self.pupsub
            elif self.symtype==1:
                useysubap=(useysubap+1)/2
                pupfn=self.pupsub[:useysubap]
            else:
                usexsubap=(usexsubap+1)/2
                useysubap=(useysubap+1)/2
                pupfn=self.pupsub[:useysubap,:usexsubap]
            for i in range(useysubap):
                for j in range(usexsubap):
                    for k in range(self.phasesize):
                        for l in range(self.phasesize):
                            indx=l+self.phasesize*k+self.phasesize**2*(j+usexsubap*i)
                            if pupfn[i,j,k,l]==1 and self.subflag[i,j]==1:
                                self.puparr[indx/8]=self.puparr[indx/8] | (1<<(indx%8))
        #finished getting it in bit form... so...
        #now check that the pupil has been created in a legal symmetrical form - if not, warn user.
        if self.symtype==1:
            tmppupsub=self.pupsub.copy()
            tmppupsub[(self.nsubx+1)/2:,]=self.pupsub[self.nsubx/2-1::-1]
            if numpy.sum((1-(tmppupsub==self.pupsub)).flat)>0:
                print "WARNING: pupil mask is not subap-symmetric about the x axis - the FPGA will be using a slightly different pupil mask"
        elif self.symtype==2:
            tmppupsub=self.pupsub.copy()
            tmppupsub[(self.nsubx+1)/2:,]=self.pupsub[self.nsubx/2-1::-1]
            tmppupsub[:,(self.nsubx+1)/2:,]=tmppupsub[:,self.nsubx/2-1::-1]
            if numpy.sum((1-(tmppupsub==self.pupsub)).flat)>0:
                print "WARNING: pupil mask is not subap-symmetric about the x axis - the FPGA will be using a slightly different pupil mask"

        tmparr=util.arrayFromArray.arrayFromArray(self.fpgaarr,(4*1024*1024,),numpy.uint8)
        savedarr=tmparr.copy()
        tmparr[:,]=self.puparr[:,]
        #now load the data to the fpga.
        print "Loading pupil function into FPGA..."
        fpid=self.fpid
        fpga.writeReg(fpid,0x2,6)#stop pipe
        fpga.writeAddr(fpid,tmparr,1)#input address
        fpga.writeReg(fpid,4*1024*1024/8,2)#size (quad words)
        fpga.writeReg(fpid,0,4)#size to write (0 bytes).
        fpga.writeReg(fpid,64,6)#reset input/output fifos.
        fpga.writeReg(fpid,8,512)#set pipe to write pupil fn.
        addr=fpga.readReg(fpid,525)
        print "Loading Pupil: Current QDR address is: %s, setting to zero"%str(addr)
        fpga.writeReg(fpid,0,525)
        fpga.writeReg(fpid,1,6)#set it going.
        while 1:#should only print a few messages here at most before done.
            time.sleep(0.008)
            addr=fpga.readReg(fpid,525)
            print "Loading Pupil: Writing to address %s"%str(addr)
            if addr==0 or addr==2**19:
                break
        print "FPGA QDR memory filled with pupil function"
        tmparr[:,]=savedarr#copy data back.
        time.sleep(0.01)
        fpga.writeReg(fpid,0,512)#unset the QDR write pipe.

    def loadSeed(self):
        """Load the random number generator seed (for readout noise).  This is taken from the Cray mta_test.c example, converted to python.  Note, this overwrites any data in reorderedPhs..."""
        print "Loading FPGA seed"
        defaultSeed=4357L
        int32_mask=0xffffffff
        multiplier=1812433253L #Don Knuth, Vol 2
        seedarr=util.arrayFromArray.arrayFromArray(self.fpgaarr,(624,),numpy.int32)#really should be UInt32, but haven't recompiled fpgamodule.c to cope yet.  Seed different, but so what!.
        savedarr=seedarr.copy()
        s=defaultSeed
        seedarr[0]=s&int32_mask
        for i in range(1,624):
            tmp=(multiplier*(seedarr[i-1]^(seedarr[i-1]>>30))+i)
            seedarr[i]=tmp&int32_mask
        #now load the data to the fpga.
        print "Random seed memory array created, loading into FPGA..."
        fpid=self.fpid
        fpga.writeReg(fpid,0x2,6)#stop pipe
        fpga.writeAddr(fpid,seedarr,1)#input address
        fpga.writeReg(fpid,624/2,2)#size (quad words)
        fpga.writeReg(fpid,0,4)#size to write (0 bytes).
        fpga.writeReg(fpid,64,6)#reset input/output fifos.
        fpga.writeReg(fpid,4,512)#set pipe to write QDR.
        fpga.writeReg(fpid,1,6)#set it going.
        time.sleep(0.01)
        print "Written seed to FPGA"
        seedarr[:,]=savedarr
        fpga.writeReg(fpid,0,512)#unset the QDR write pipe.
        
    def loadPoisson(self,calsource):
        """Load the QDR memory with appropriate random variables.
        For <2ppp, have 256 bins, each with 512 byte entries. (64 qwords)
        For <8ppp, have 384 bins, each with 1024 byte entries. (128 qwords)
        For <32ppp, have 768 bins, each with 2048 byte entries. (256 qwords)
        For >=32ppp, have gaussian. (262144 qwords, 2MB, 1048576 short entries)
          --addressing scheme: 20 bits (only 19 valid for QDR).
          --01234567890123456789
          --000000ppppppppcccccc     <2 (14 bits)
          --0000aapppppppccccccc    2-8 (16 bits)
          --00aappppppppcccccccc   8-32 (18 bits)
          --01iiiiiiiiiiiiiiiiii   Gaussian.
          --(aa here means that at least one of these numbers will be 1.
          --p represents bits coming from the light level, c represents bits
          --coming from the bin, ie the dpbm output)
        Note, this overwrites any data in the reorderedPhs array...
        """
        print "Loading FPGA Poisson seed"
        qdrmem=util.arrayFromArray.arrayFromArray(self.fpgaarr,(4*1024*1024,),numpy.uint8)#everything
        qdrmemp=util.arrayFromArray.arrayFromArray(self.fpgaarr,(2*1024*1024,),numpy.uint8)#poisson part
        qdrmems=util.arrayFromArray.arrayFromArray(self.fpgaarr[2*1024*1024:,],(1024*1024,),numpy.int16)#gaussian part - signed 8.8 format.
        savedarr=qdrmem.copy()
        for i in range(256):#less than 2ppp, fill the array.
            ppp=i*2**-7#the mean light level of this bin.
            #RandomArray.poisson(ppp,(512,)).astype(numpy.UInt8)
            if calsource==0:
                qdrmemp[i*512:i*512+512]=util.poisson.mypoissondist(ppp,512)
            else:
                qdrmemp[i*512:i*512+512]=int(ppp)
        for i in range(384):
            #2-8ppp.
            ppp=2+i*2**-6
            if calsource==0:
                qdrmemp[256*512+i*1024:256*512+i*1024+1024]=util.poisson.mypoissondist(ppp,1024)
            else:
                qdrmemp[256*512+i*1024:256*512+i*1024+1024]=int(ppp)
        for i in range(768):
            ppp=8+i*2**-5
            if calsource==0:
                qdrmemp[256*512+384*1024+i*2048:256*512+384*1024+i*2048+2048]=util.poisson.mypoissondist(ppp,2048)
            else:
                qdrmemp[256*512+384*1024+i*2048:256*512+384*1024+i*2048+2048]=int(ppp)
        if calsource==0:
            qdrmems[:,]=(numpy.random.standard_normal(1024*1024)*2**8).astype(numpy.int16)
        else:
            qdrmems[:,]=0
        self.savedQDRMem=qdrmem.copy()

        #now we have the memory done, should load it into the FPGA.
        print "QDR memory array created, loading into FPGA..."
        fpid=self.fpid
        fpga.writeReg(fpid,0x2,6)#stop pipe
        fpga.writeAddr(fpid,qdrmem,1)#input address
        fpga.writeReg(fpid,4*1024*1024/8,2)#size (quad words)
        fpga.writeReg(fpid,0,4)#size to write (0 bytes).
        fpga.writeReg(fpid,64,6)#reset input/output fifos.
        fpga.writeReg(fpid,1,512)#set pipe to write QDR.
        addr=fpga.readReg(fpid,525)
        print "Loading poisson RV to FPGA: Current QDR address is %s, setting to zero"%str(addr)
        fpga.writeReg(fpid,0,525)
        fpga.writeReg(fpid,1,6)#set it going.
        while 1:#should only print a few messages here at most before done.
            time.sleep(0.008)
            addr=fpga.readReg(fpid,525)
            print "Loading Poisson RV: Writing to address %s"%str(addr)
            if addr==0 or addr==2**19:
                break
        print "QDR memory filled"
        qdrmem[:,]=savedarr
        fpga.writeReg(fpid,0,512)#unset the QDR write pipe.

    def loadBinCtrl(self):
        """Load control vectors for binning into the FPGA"""
        #now set up the binning control...
        tmp=numpy.zeros((32,),numpy.int8)
        cnt=0
        #first do x...
        for i in range(32):
            cnt+=self.nimg
            if cnt>=self.fftsize:
                excess=cnt-self.fftsize
                use=self.nimg-excess
                tmp[i]=(use&0x3f)|0x40
                cnt=excess
            else:
                tmp[i]=(self.nimg&0x3f)
        larr=util.arrayFromArray.arrayFromArray(tmp,(4,),numpy.int64)
        for i in range(4):
            fpga.writeReg(self.fpid,larr[i],536+i)
        cnt=0
        #now for the y...
        for i in range(32):
            cnt+=self.nimg
            if cnt>=self.fftsize:
                excess=cnt-self.fftsize
                use=self.nimg-excess
                tmp[i]=(use&0x3f)|0x40
                cnt=excess
            else:
                tmp[i]=(self.nimg&0x3f)
        for i in range(4):
            fpga.writeReg(self.fpid,larr[i],540+i)
        #have now finished setting up the binning control.
    def loadCenMask(self):
        """Load the centroid mask into the FPGA registers..."""
        iarr=numpy.zeros((32,),numpy.uint32)
        larr=util.arrayFromArray.arrayFromArray(iarr,(16,),numpy.int64)
        #Largest cent mask to be loaded is 32x32 pixels.  So, have to map
        #out cent mask onto an array this size.
        pos=0
        cenmask=self.cenmask.astype(numpy.int64)
        for y in range(self.nimg):
            for x in range(self.nimg):
                iarr[pos/32]=iarr[pos/32] | (cenmask[y,x]<<(pos%32))
                pos+=1
        #print larr
        for i in range(16):
            fpga.writeReg(self.fpid,larr[i],544+i)
        
    def loadFPGARegAddr(self):
        """Initialise the FPGA registers"""
        print "Loading FPGA address registers"
        fpid=self.fpid
        fpga.writeReg(fpid,0x2,6)#stop pipe
        if self.doPartialFPGA:
            fpga.writeAddr(fpid,self.fpgaInArr,1)#input address
            fpga.writeReg(fpid,self.fittedSubaps*self.nIntegrations*self.phasesize*self.phasesize*4/8,2)#size (quad words)
            fpga.writeAddr(fpid,self.fpgaOutArr,3)#output address
            fpga.writeReg(fpid,self.fittedSubaps,4)#size to write (in bytes).
        else:
            fpga.writeAddr(fpid,self.reorderedPhs,1)#input address
            fpga.writeReg(fpid,self.nsubx*self.nsubx*self.nIntegrations*self.phasesize*self.phasesize*4/8,2)#size (quad words)
            fpga.writeAddr(fpid,self.outputData,3)#output address
            fpga.writeReg(fpid,self.nsubx*self.nsubx,4)#size to write (in bytes).
    def loadFPGADimData(self):
        """Break loading of registers up into a few sections..."""
        print "Loading FPGA dimensions"
        fpid=self.fpid
        fpga.writeReg(fpid,self.nIntegrations,519)#number of integrations
        fpga.writeReg(fpid,self.nimg,520)#centroid subap x size
        fpga.writeReg(fpid,self.nimg,521)#centroid subap y size
        fpga.writeReg(fpid,self.fftsize,522)#fft size
        fpga.writeReg(fpid,self.phasesize,523)#phase size
        fpga.writeReg(fpid,2**15/self.fftsize,526)#inv fft size
        #fpga.writeReg(fpid,self.subapsize,527)#scale/norm subapsize
        fpga.writeReg(fpid,self.nimg,532)#centroid subap x size (again)
        fpga.writeReg(fpid,self.nimg,533)#centroid subap y size (again)
        fpga.writeReg(fpid,(self.nsubx&0x3ff)|((self.nsubx&0x3ff)<<10),534)#write the number of subaps
        fpga.writeReg(fpid,self.symtype,535)#write the symmetry type
    def loadFPGASourceData(self,calsource):
        print "Loading FPGA source info"
        fpid=self.fpid
        fpga.writeReg(fpid,self.sigFPGA,524)#user scale
        fpga.writeReg(fpid,(1-calsource)*self.fpgaSkybrightness,528)#sky brightness in 16.10 fixed point format per pupil pixel in each subap...
        fpga.writeReg(fpid,(1-calsource)*int(self.noiseFloor),529)#threshold value to remove
        fpga.writeReg(fpid,(1-calsource)*self.fpga_readbg,530)#mean readout noise (const level added to signal)
        fpga.writeReg(fpid,(1-calsource)*self.fpga_readnoise,531)#standard deviation of readout noise (RMS readout noise)

    def tilt(self,x,y):
        return y

    def computeHLL(self,phase):
        """Compute high light level image from phase, with nsubx subapertures
        in x and y directions.
        oversample is the factor used for oversampling the FFT.
        convarr is an optional array image (2 or 4d) which gets convolved
        with the SHS spots - eg an airy disc or LGS spot elongations...
        If an airy disc, create it with xoff=yoff=0.5.
        """
        nsubx=self.nsubx
        pupfn=self.pupfn
        oversample=self.oversamplefactor
        self.wfsn=n=phase.shape[0]/nsubx
        nfft=n*oversample
        subimg=numpy.zeros((nsubx,nsubx,nfft,nfft),"d")
        #fftwPlan=cmod.mkimg.setup(subimg[0,0])
        pupsub=numpy.ones((nsubx,nsubx,n,n),"d")
        if type(pupfn)!=type(None):
            if pupfn.shape!=phase.shape:
                print "ERROR: centroid.py - pupil function shape not equal to phase shape."
            for i in range(nsubx):
                for j in range(nsubx):
                    pupsub[i,j]=pupfn[i*n:(i+1)*n,j*n:(j+1)*n].astype("d")
        phase=numpy.array(phase)
        tmp=0.5*n/nfft*2*numpy.pi
        xtiltfn=((numpy.fromfunction(self.tilt,(n,n))-float(n)/2.+0.5)/float(n)).astype("d")# subap tilt fn

        tiltfn=tmp*(xtiltfn+numpy.transpose(xtiltfn))
        for i in range(nsubx):
            for j in range(nsubx):
                phs=phase[i*n:i*n+n,j*n:j*n+n].astype("d")
                phs-=tiltfn
                #cmod.mkimg.mkimg(fftwPlan,phs,subimg[i,j],pupsub[i,j])
                self.makeImage(phs,subimg[i,j],pupsub[i,j])
                subimg[i,j]=util.flip.fliparray(subimg[i,j])
                if type(self.psf)!=type(None):
                    if len(self.psf.shape)==4:
                        ca=self.psf[i,j]
                    else:
                        ca=self.psf
                    subimg[i,j]=self.conv(subimg[i,j],ca)
        self.subimg=subimg
        return subimg#long exp image.

    def conv(self,img1,img2):
        """Convolve images, eg subap image with a PSF
        e.g. LGS spot shape or airy disc.
        Note, if using an airy disc form img2, create it with xoff=yoff=0.5"""
##         nfft=self.wfs_nfft
##         temp=numpy.zeros((nfft,nfft),numpy.float64)
##         temp[nfft/2-1,nfft/2-1]=1.
##         temp[nfft/2,nfft/2]=1.
##         temp[nfft/2-1,nfft/2]=1.
##         temp[nfft/2,nfft/2-1]=1.
##         wfsimg=fliparray(temp)
        #import FFT
        temp1=numpy.fft.rfft2(img1)
        temp2=numpy.fft.rfft2(util.flip.fliparray2(img2))
        convimg=numpy.fft.irfft2(temp1*temp2)
        return convimg

        
        

    def makeshimg(self,subimg):
        shape=subimg.shape
        self.tilenoiseless=numpy.zeros((shape[0]*shape[2],shape[1]*shape[3]),subimg.dtype)
        for i in range(shape[0]):
            for j in range(shape[1]):
                self.tilenoiseless[i*shape[2]:(i+1)*shape[2],j*shape[3]:(j+1)*shape[3]]=subimg[i,j]
        return self.tilenoiseless
    
    def calcCents(self,subimg,ninteg=1):
        binfactor=self.binfactor
        nsubx=subimg.shape[0]
        nfft=subimg.shape[2]
        nimg=nfft/binfactor
        self.nimg=nimg
        n=nfft/self.oversamplefactor
        nsubs=nsubx*nsubx
        bimg=numpy.zeros((nsubx,nsubx,nimg,nimg),"d")
        self.bimg=bimg
        cent=numpy.zeros((nsubs*2,),"d")
        self.centx=util.arrayFromArray.arrayFromArray(cent[:nsubs],(nsubx,nsubx),cent.dtype)
        self.centy=util.arrayFromArray.arrayFromArray(cent[nsubs:,],(nsubx,nsubx),cent.dtype)
        self.tile=numpy.zeros((nimg*nsubx,nimg*nsubx),subimg.dtype)
        self.photPerSubap=numpy.zeros((nsubx,nsubx),"d")
        indx = (numpy.array(range(nimg),numpy.float64))-float(nimg/2)+0.5
        for i in range(nsubx):
            for j in range(nsubx):
                cmod.binimg.binimg(subimg[i,j],bimg[i,j])
                totsig=numpy.sum(numpy.sum(bimg[i,j]))
                if type(self.subarea)==type(None):
                    nphs=1.
                else:
                    nphs=float(self.subarea[i,j])/(n*n)#fraction of active pxls.
                if totsig>0:
                    if type(self.sig)==type(0.0):
                        bimg[i,j]*=self.sig*nphs/totsig
                    else:
                        bimg[i,j]*=self.sig[i,j]/totsig
                bimg[i,j]+=self.skybrightness*nphs
                for k in range(ninteg):#perform several integrations...
                    tmpimg=bimg[i,j].copy()
                    if self.addPoisson:
                        cmod.imgnoise.shot(tmpimg,tmpimg)
                    if self.readnoise>0 or self.readbg>0:
                        tmpimg+=numpy.random.normal(self.readbg,self.readnoise,tmpimg.shape)
                    tmpimg[:,]=numpy.where(tmpimg<self.noiseFloor,0,tmpimg-self.noiseFloor)
                bimg[i,j]=tmpimg
                self.tile[i*nimg:(i+1)*nimg,j*nimg:(j+1)*nimg]=bimg[i,j]
                totsig=numpy.sum(numpy.sum(bimg[i,j]))
                self.photPerSubap[i,j]=totsig
                if totsig==0:
                    totsig=1.
                if self.centroidPower==None:
                    self.centx[i,j]=numpy.sum(numpy.sum(bimg[i,j],0)*indx)/totsig
                    self.centy[i,j]=numpy.sum(numpy.sum(numpy.transpose(bimg[i,j]),0)*indx)/totsig
                else:
                    self.centx[i,j]=numpy.sum(numpy.sum(bimg[i,j]**self.centroidPower,0)*indx)/totsig
                    self.centy[i,j]=numpy.sum(numpy.sum(numpy.transpose(bimg[i,j]**self.centroidPower),0)*indx)/totsig

        self.cent=cent
        if self.warnOverflow!=None:
            if max(self.tile.flat)>self.warnOverflow:
                print "Max value in CCD is",max(self.tile.flat)
        if self.printmax:
            print "Max value in CCD is",max(self.tile.flat)
            print "Mean signal in CCD is",numpy.average(self.tile.flat)
            
        return self.cent

    def run(self,phase,ninteg=1):
        hll=self.computeHLL(phase)
        simg=self.makeshimg(hll)
        c=self.calcCents(hll,ninteg)
        return self

    def calc(self,phase):
        return self.calcCents(self.computeHLL(phase))

    def getCentral(self,npxls=4,offset=0):
        """Return an image of SH spots, but only the central few spots are shown.  This can be useful if you have a system with a large number of pxls per subap, and a fairly flat wavelength, and wish to view eg in gist..."""
        nsubx=self.nsubx
        nimg=self.nimg
        s=nimg/2-npxls/2+offset
        e=s+npxls
        img=numpy.zeros((npxls*nsubx,npxls*nsubx),"d")
        for i in range(nsubx):
            for j in range(nsubx):
                img[i*npxls:(i+1)*npxls,j*npxls:(j+1)*npxls]=self.bimg[i,j,s:e,s:e]
        return img
    def computePxlToRad(self,n):
        c=centroid(1,oversamplefactor=self.oversamplefactor,binfactor=self.binfactor)
        phs=numpy.zeros((n,n),"d")
        phs[:,]=numpy.arange(n).astype("d")/(n-1)
        cc=c.calc(phs)
        self.convFactor=abs(1./cc[0])#multiply the centroids (in pixels) by this value to get the mean phase slopes (centroids) in radians.


    def computeNoiseCovariance(self,niters=10,wfsn=8,convertToRad=0):
        """puts flat phase in, and computes centroids for several iterations,
        and then computes covariance."""
        if type(self.pupfn)!=type(None):
            phase=numpy.zeros(self.pupfn.shape,numpy.float64)
        else:
            phase=numpy.zeros((self.nsubx*wfsn,self.nsubx*wfsn),numpy.float64)
        n=self.nsubx*self.nsubx*2
        c2=numpy.zeros((n,),numpy.float64)
        csum=numpy.zeros((n,),numpy.float64)
        for i in range(niters):
            self.run(phase)
            csum+=self.cent
            c2+=self.cent*self.cent
        mean=csum/n
        var=c2/n-mean*mean
        convFactor=1.
        if convertToRad:
            self.computePxlToRad(wfsn)
            convFactor=self.convFactor**2
##             #this is a bit of a wasteful way of doing it but never mind!
##             #First, compute a centroid with a slope of xxx radians.
##             #This will then tell us that a pxl of yyy value equals
##             #xxx radians.
##             #Do this using noiseless centroider.
##             rn=self.readnoise
##             ap=self.addPoisson
##             rb=self.readbg
##             nf=self.noiseFloor
##             sb=self.skybrightness
##             sa=self.subarea
##             pf=self.pupfn
##             nsx=self.nsubx
##             self.readnoise=0.
##             self.addPoisson=0
##             self.readbg=0.
##             self.noiseFloor=0.
##             self.skybrightness=0.
##             self.subarea=None
##             self.pupfn=None
##             phase[:,]=numpy.arange(phase.shape[1]).astype("d")/(phase.shape[1]-1)*self.nsubx#range from 0-1 radian of phase across each subaperture.
##             cents=self.calc(phase)#should all be same?
##             #So a 1 radian phase tilt across a subap gives a centroid value of cent[0].
##             #So an RMS centroid value of XXX corresponds to a RMS phase value of XXX/cent[0].  
##             #So, Var of XXX corresponds to value in rad^2 of XXX/cent[0]**2
##             self.convFactor=1./cents[0]/cents[0]
##             #now reinstate the noise...
##             self.readnoise=rn
##             self.addPoisson=ap
##             self.readbg=rb
##             self.noiseFloor=nf
##             self.skybrightness=sb
##             self.subarea=sa
##             self.pupfn=pf
            
        return var*convFactor
        
    def magicShackHartmann(self):
        """Taken from TB - just takes slope measurement."""
        nsubx=self.nsubx
        pupfn=self.pupfn
        n=self.pupfn.shape[0]/nsubx
        if self.magicSlopes==None:
            self.magicSlopes=self.magicSHSlopes()
        magicSlopes=self.magicSlopes

        #now need to reorder... first average the different integration times
        #Actually, maybe we should just use the latest, since we're magic... oh well.
        phs=(self.reorderedPhs.sum(2)/self.reorderedPhs.shape[2])[:,:,:,:self.phasesize]
        #phs.shape=nsubx,nsubx,phasesize,phasesize_v
        

        #tmp1=phs*magicSlopes#*R0_LAMBDA/WFS_LAMBDA
        #tmp2=phs*magicSlopes.transpose()#*R0_LAMBDA/WFS_LAMBDA
        cent=self.outputData#numpy.zeros((nsubx**2,2),numpy.float32)
        for i in range(nsubx):
            for j in range(nsubx):
                if self.subflag[i,j]:
                    cent[i,j,0]=(phs[i,j]*magicSlopes[i*n:(i+1)*n,j*n:(j+1)*n]).sum()
                    cent[i,j,1]=(phs[i,j]*magicSlopes.transpose()[i*n:(i+1)*n,j*n:(j+1)*n]).sum()
        return cent

    def magicSHSlopes(self):
        npup=self.pupfn.shape[0]
        slop=numpy.zeros((npup,npup),numpy.float32)
        slop[:]=numpy.arange(npup)+1.
        slop*=self.pupfn
        n=self.pupfn.shape[0]/self.nsubx
        #self.subflag=numpy.zeros((self.nsubx,self.nsubx),numpy.int32)
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                #if self.pupfn[i*n:(i+1)*n,j*n:(j+1)*n].sum()>=self.wfs_minarea*n*n:
                #   self.subflag[i,j]=1
                if self.subflag[i,j]>0:
                    p=slop[i*n:(i+1)*n,j*n:(j+1)*n].sum()/float(self.pupfn[i*n:(i+1)*n,j*n:(j+1)*n].sum())
                    slop[i*n:(i+1)*n,j*n:(j+1)*n]-=p*self.pupfn[i*n:(i+1)*n,j*n:(j+1)*n]
                    s2=(slop[i*n:(i+1)*n,j*n:(j+1)*n]**2).sum()
                    slop[i*n:(i+1)*n,j*n:(j+1)*n]/=-numpy.sqrt(s2)
                else:
                    slop[i*n:(i+1)*n,j*n:(j+1)*n]*=0.
        return slop

    def calibrateSHSUnique(self,control={"cal_source":1,"useFPGA":0,"useCell":0,"useCmod":1}):
        if self.linearSteps==None:
            return
        print "Calibrating centroids (all subaps treated differently)"
        steps=self.linearSteps
        self.linearSteps=None
        #create a (fairly large) array to store the data.
        self.calibrateData=numpy.zeros((2,self.nsubx,self.nsubx,steps),numpy.float32)
        #compute the phase slopes to use...
        stepList=(numpy.arange(steps)-steps/2.+0.5)/(steps/2.-0.5)*numpy.pi*self.stepRangeFrac

        self.calibrateSteps=stepList.astype(numpy.float32)
        c=control["cal_source"]
        #control["cal_source"]=1
        control["cal_source"]=0
        if self.centcmod!=None:
            self.centcmod.update(util.centcmod.ADDPOISSON,0)
            self.centcmod.update(util.centcmod.READNOISE,0)
        else:
            raise Exception("Not yet sorted for non-cmod type things")
        for i in range(steps):
            #put a slope into it
            tilt=-numpy.arange(self.phasesize)*stepList[i]
            self.reorderedPhs[:,:,:]=tilt
            #compute centroids (x)
            self.runCalc(control)
            #Now take the x centroids and store...
            self.calibrateData[0,:,:,i]=self.centx

            #And now repeat for the y centroids
            self.reorderedPhs[:,:,:]=tilt[None].transpose()
            self.runCalc(control)
            self.calibrateData[1,:,:,i]=self.centy
        if self.centcmod!=None:#restore the original values.
            self.centcmod.update(util.centcmod.ADDPOISSON,self.centcmod.addPoisson)
            self.centcmod.update(util.centcmod.READNOISE,self.centcmod.readnoise)
        control["cal_source"]=c
        self.linearSteps=steps
        #self.calibrateDataOrig=self.calibrateData.copy()
        #Now compute the bounds for which this is single valued (since, when cent gets too close to edge, it starts wrapping round).
        self.calibrateBounds=numpy.zeros((2,self.nsubx,self.nsubx,2),numpy.int32)
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                #x
                self.calibrateBounds[0,i,j,0]=numpy.argmin(self.calibrateData[0,i,j])
                self.calibrateBounds[0,i,j,1]=numpy.argmax(self.calibrateData[0,i,j])
                #and y
                self.calibrateBounds[1,i,j,0]=numpy.argmin(self.calibrateData[1,i,j])
                self.calibrateBounds[1,i,j,1]=numpy.argmax(self.calibrateData[1,i,j])
        #Also need to check that the calibrateData[:,:,:,i] is increasing always - otherwise the interpolation won't work.  What should we do if its not increasing???
        if self.calNCoeff==0:
            for i in range(self.nsubx):
                for j in range(self.nsubx):
                    for k in range(self.calibrateBounds[0,i,j,0],self.calibrateBounds[0,i,j,1]):
                        if self.calibrateData[0,i,j,k]>self.calibrateData[0,i,j,k+1]:
                            val=(self.calibrateData[0,i,j,k-1]+self.calibrateData[0,i,j,k+1])/2.
                            print "Forcing SHS calibration for point (%d,%d) step %d from %g to %g"%(i,j,k,self.calibrateData[0,i,j,k],val)
                            self.calibrateData[0,i,j,k]=val
                    for k in range(self.calibrateBounds[1,i,j,0],self.calibrateBounds[1,i,j,1]):
                        if self.calibrateData[1,i,j,k]>self.calibrateData[1,i,j,k+1]:
                            val=(self.calibrateData[1,i,j,k-1]+self.calibrateData[1,i,j,k+1])/2.
                            print "Forcing SHS calibration for point (%d,%d) step %d from %g to %g"%(i,j,k,self.calibrateData[1,i,j,k],val)
                            self.calibrateData[1,i,j,k]=val
        print "Finished calibrating centroids"

    def applyCalibrationUnique(self,data=None):
        """Uses the calibration, to replace data with a calibrated version of data.
        Data.shape should be nsubx,nsubx,2
        Typically it will be self.outputData
        """
        if data==None:
            data=self.outputData
        if self.calNCoeff==0:
            for i in range(self.nsubx):
                for j in range(self.nsubx):
                    if self.subflag[i,j]:
                        cx,cy=data[i,j]#the x,y centroids.
                        if cx>self.calibrateData[0,i,j,self.calibrateBounds[0,i,j,1]] or cx<self.calibrateData[0,i,j,self.calibrateBounds[0,i,j,0]]:
                            print "Warning: x centroid at %d,%d with value %g is outside the calibrated bounds"%(i,j,cx)
                        if cy>self.calibrateData[1,i,j,self.calibrateBounds[1,i,j,1]] or cy<self.calibrateData[1,i,j,self.calibrateBounds[1,i,j,0]]:
                            print "Warning: y centroid at %d,%d with value %g is outside the calibrated bounds"%(i,j,cy)
                        data[i,j,0]=numpy.interp([cx],self.calibrateData[0,i,j,self.calibrateBounds[0,i,j,0]:self.calibrateBounds[0,i,j,1]+1],self.calibrateSteps[self.calibrateBounds[0,i,j,0]:self.calibrateBounds[0,i,j,1]+1])[0]
                        #print "applyCalibration error %d %d 0 %g %d %d"%(i,j,cx,self.calibrateBounds[0,i,j,0],self.calibrateBounds[0,i,j,1]+1)
                        data[i,j,1]=numpy.interp([cy],self.calibrateData[1,i,j,self.calibrateBounds[1,i,j,0]:self.calibrateBounds[1,i,j,1]+1],self.calibrateSteps[self.calibrateBounds[1,i,j,0]:self.calibrateBounds[1,i,j,1]+1])[0]
                        #print "applyCalibration error %d %d 1 %g %d %d"%(i,j,cy,self.calibrateBounds[1,i,j,0],self.calibrateBounds[1,i,j,1]+1)
        else:#calibration using interpolation
            klist=range(self.calNCoeff)
            for i in range(self.nsubx):
                for j in range(self.nsubx):
                    if self.subflag[i,j]:
                        resx=0.
                        resy=0.
                        xc=1.
                        yc=1.
                        cx,cy=data[i,j]
                        for k in klist:
                            resx+=xc*self.calCoeff[i,j,0,k]
                            xc*=cx
                            resy+=yc*self.calCoeff[i,j,1,k]
                            yc*=cy
                        data[i,j]=resx,resy

    def calibrateSHSIdentical(self,control={"cal_source":1,"useFPGA":0,"useCell":0,"useCmod":1}):
        if self.linearSteps==None:
            return
        print "Calibrating centroids (identical subap pupil functions treated same)"
        steps=self.linearSteps
        self.linearSteps=None
        #create a (fairly large) array to store the data.
        self.calibrateData=numpy.zeros((2,self.nsubx,self.nsubx,steps),numpy.float32)
        self.calMaskType=numpy.zeros((self.nsubx,self.nsubx),numpy.int32)
        self.calDataDict={}
        typ=1
        d={0:[]}
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                if self.subflag[i,j]:
                    keys=d.keys()
                    typ=None
                    for key in keys:
                        if numpy.alltrue(d[key]==self.pupfn[i*self.phasesize:(i+1)*self.phasesize,j*self.phasesize:(j+1)*self.phasesize]):
                            typ=key
                            break
                    if typ==None:#new type needed
                        typ=max(keys)+1
                        d[typ]=self.pupfn[i*self.phasesize:(i+1)*self.phasesize,j*self.phasesize:(j+1)*self.phasesize]
                    self.calMaskType[i,j]=typ


        #compute the phase slopes to use...
        stepList=(numpy.arange(steps)-steps/2.+0.5)/(steps/2.-0.5)*numpy.pi*self.stepRangeFrac
        self.calibrateSteps=stepList.astype(numpy.float32)
        c=control["cal_source"]
        #control["cal_source"]=1
        control["cal_source"]=0
        if self.centcmod!=None:
            self.centcmod.update(util.centcmod.ADDPOISSON,0)
            self.centcmod.update(util.centcmod.READNOISE,0)
        else:
            raise Exception("Not yet sorted for non-cmod type things")
        #TODO: the best way of doing this is to use cal_source==0, and set addPoisson=0, readnoise=0, and that way the backgrounds are correct.
        for i in range(steps):
            #put a slope into it
            tilt=-numpy.arange(self.phasesize)*stepList[i]
            self.reorderedPhs[:,:,:]=tilt
            #compute centroids (x)
            self.runCalc(control)
            #Now take the x centroids and store...
            self.calibrateData[0,:,:,i]=self.centx

            #And now repeat for the y centroids
            self.reorderedPhs[:,:,:]=tilt[None].transpose()
            self.runCalc(control)
            self.calibrateData[1,:,:,i]=self.centy
        if self.centcmod!=None:#restore the original values.
            self.centcmod.update(util.centcmod.ADDPOISSON,self.centcmod.addPoisson)
            self.centcmod.update(util.centcmod.READNOISE,self.centcmod.readnoise)
        control["cal_source"]=c
        self.linearSteps=steps
        #self.calibrateDataOrig=self.calibrateData.copy()
        #Now compute the bounds for which this is single valued (since, when cent gets too close to edge, it starts wrapping round).
        self.calibrateBounds=numpy.zeros((2,self.nsubx,self.nsubx,2),numpy.int32)
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                if self.subflag[i,j]:
                    if not self.calDataDict.has_key(self.calMaskType[i,j]):
                        # Create a new store for this type of subap
                        xf=numpy.argmin(self.calibrateData[0,i,j])
                        xt=numpy.argmax(self.calibrateData[0,i,j])+1
                        yf=numpy.argmin(self.calibrateData[1,i,j])
                        yt=numpy.argmax(self.calibrateData[1,i,j])+1
                        xc=self.calibrateData[0,i,j,xf:xt].copy()
                        yc=self.calibrateData[1,i,j,yf:yt].copy()
                        xr=self.calibrateSteps[xf:xt]
                        yr=self.calibrateSteps[yf:yt]
                        self.calDataDict[self.calMaskType[i,j]]=CalData(xc,xr,yc,yr)
                        cd=self.calDataDict[self.calMaskType[i,j]]
                        #mod=0
                        for k in range(cd.xc.shape[0]-1):
                            if cd.xc[k]>=cd.xc[k+1]:
                                #mod=1
                                #print cd.xc
                                #import util.FITS
                                #util.FITS.Write(cd.xc,"tmp.fits")
                                if k==0:
                                    val=cd.xc[k+1]-(cd.xc[k+2]-cd.xc[k+1])
                                else:
                                    val=(cd.xc[k-1]+cd.xc[k+1])/2.
#                                 if val==cd.xc[k+1]:
#                                     if val>0:
#                                         val*=0.99999
#                                     else:
#                                         val*=1.00001
                                print "Forcing SHS x calibration for point (%d,%d) step %g from %g to %g"%(i,j,cd.xr[k],cd.xc[k],val)
                                cd.xc[k]=val
                        notsame=numpy.nonzero(cd.xc[1:]!=cd.xc[:-1])[0]
                        cd.xc=cd.xc[notsame]
                        cd.xr=cd.xr[notsame]
                        for k in range(cd.yc.shape[0]-1):
                            if cd.yc[k]>=cd.yc[k+1]:
                                if k==0:
                                    val=cd.yc[k+1]-(cd.yc[k+2]-cd.yc[k+1])
                                else:
                                    val=(cd.yc[k-1]+cd.yc[k+1])/2.
#                                 if val==cd.yc[k+1]:
#                                     if val>0:
#                                         val*=0.99999
#                                     else:
#                                         val*=1.00001
                                print "Forcing SHS y calibration for point (%d,%d) step %g from %g to %g"%(i,j,cd.yr[k],cd.yc[k],val)
                                cd.yc[k]=val
                        notsame=numpy.nonzero(cd.yc[1:]!=cd.yc[:-1])[0]
                        cd.yc=cd.yc[notsame]
                        cd.yr=cd.yr[notsame]
                        #if mod:
                        #    util.FITS.Write(cd.xc,"tmp.fits")
                    # and store this subap index.
                    cd=self.calDataDict[self.calMaskType[i,j]]
                    cd.indx.append(i*self.nsubx+j)
        for k in self.calDataDict.keys():
            cd=self.calDataDict[k]
            cd.xindx=numpy.array(cd.indx).astype(numpy.int32)*2
            cd.yindx=cd.xindx+1
        print "Finished calibrating centroids"

    def applyCalibrationIdentical(self,data=None):
        """Uses the calibration, to replace data with a calibrated version of data.
        Data.shape should be nsubx,nsubx,2
        Typically it will be self.outputData
        """
        if data==None:
            data=self.outputData
        if numpy.any(numpy.isnan(data)):
            print "Warning -nan prior to applyCalibrationIdentical"
        rdata=data.ravel()
        for k in self.calDataDict.keys():
            cd=self.calDataDict[k]
            x=numpy.take(rdata,cd.xindx)
            y=numpy.take(rdata,cd.yindx)
            #print "x,y",x,y
            if numpy.any(x>cd.xc[-1]) or numpy.any(x<cd.xc[0]):
                print "Warning: x centroid is outside calibrated bounds"
            if numpy.any(y>cd.yc[-1]) or numpy.any(y<cd.yc[0]):
                print "Warning: y centroid is outside calibrated bounds"
            #now put the calibrated values back (using numpy fancy indexing)
            rdata[cd.xindx]=numpy.interp(x,cd.xc,cd.xr).astype(numpy.float32)
            rdata[cd.yindx]=numpy.interp(y,cd.yc,cd.yr).astype(numpy.float32)
            indx=numpy.nonzero(numpy.isnan(rdata[cd.xindx]))[0]
            if indx.size>0:
                print "indx",indx,numpy.take(rdata[cd.xindx],indx),numpy.take(x,indx)
                print cd.xc,cd.xr
            #print "rdata",rdata[cd.xindx],rdata[cd.yindx]
        if not data.flags.c_contiguous:
            #need to copy the data back in.
            print "Warning - flattening non-contiguous centroid data"
            rdata.shape=data.shape
            data[:]=rdata
        if numpy.any(numpy.isnan(data)):
            print "Warning -nan after applyCalibrationIdentical"

    def makeCalibrationCoeffs(self):
        self.calCoeff=numpy.zeros((self.nsubx,self.nsubx,2,self.calNCoeff),numpy.float32)
        print "todo makeCalibrationCoeffs - check whether is giving best performance (calibrateSHSUnique does some adjustment to the data, so it may not) - should be use the bounds, or use the whole think unadjusted...?"
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                if self.subflag[i,j]:
                    f,t=self.calibrateBounds[0,i,j]
                    self.calCoeff[i,j,0]=numpy.polyfit(self.calibrateData[0,i,j,f:t],self.calibrateSteps[f:t],self.calNCoeff-1)[::-1]
                    self.calCoeff[i,j,1]=numpy.polyfit(self.calibrateData[1,i,j,f:t],self.calibrateSteps[f:t],self.calNCoeff-1)[::-1]


    def calibrateSHS(self,control={"cal_source":1,"useFPGA":0,"useCell":0,"useCmod":1}):
        if self.linearSteps==None:
            return
        if self.psf==None and self.correlationCentroiding==0 and self.calNCoeff==0:
            self.calibrateSHSIdentical(control)
        else:
            self.calibrateSHSUnique(control)
            if self.calNCoeff!=0:
                self.makeCalibrationCoeffs()
                if self.centcmod!=None:
                    self.centcmod.update(util.centcmod.CALCOEFF,self.calCoeff)
            else:
                if self.centcmod!=None:
                    self.centcmod.update(util.centcmod.CALDATA,(self.calibrateData,self.calibrateBounds,self.calibrateSteps))

    def applyCalibration(self,data=None):
        if self.psf==None and self.correlationCentroiding==0 and self.calNCoeff==0:
            self.applyCalibrationIdentical(data)
        else:
            if self.centcmod==None:#otherwise its been done in the c module.
                self.applyCalibrationUnique(data)


    def takeReference(self,control):
        """Measure noiseless centroid offsets for flat input.  These are then subsequently used as reference centroids.
        """
        self.reorderedPhs[:]=0
        # compute centroids (x)
        c=control["cal_source"]
        control["cal_source"]=1
        #steps=self.linearSteps
        #self.linearSteps=None
        self.runCalc(control)
        control["cal_source"]=c
        #self.linearSteps=steps
        # Now take the x centroids and store...
        self.refCents=self.outputData.copy()
        #Now, we need reference to be taken after calibration has been done.
        #So, we should only pass refs to cmod if cmod is also doing calibration (or if no calibration is being done).
        if self.centcmod!=None:
            if self.linearSteps==None or self.psf!=None or self.correlationCentroiding!=None or self.calNCoeff!=0:#no calibration done, or done in c, so ref can be done by c.
                self.centcmod.update(util.centcmod.REFCENTS,self.refCents)
                    
    def takeCorrImage(self,control):
        """If correlationCentroiding==1, but corrPattern==None, use a default SH spot pattern as the reference.
        """
        if self.correlationCentroiding:
            if self.corrPattern==None:
                c=control["cal_source"]
                control["cal_source"]=1
                steps=self.linearSteps
                self.linearSteps=None
                if self.centcmod!=None:
                    self.centcmod.update(util.centcmod.CORRELATIONCENTROIDING,0)
                else:
                    raise Exception("Not yet sorted for non-cmod type things")
                self.reorderedPhs[:]=0
                self.runCalc(control)
                self.centcmod.update(util.centcmod.CORRELATIONCENTROIDING,1)
                control["cal_source"]=c
                self.linearSteps=steps
                self.corrPatternUser=self.cmodbimg.copy()
                self.corrPatternUser/=max(self.corrPatternUser.ravel())#normalise
                self.corrPattern=util.correlation.transformPSF(self.corrPatternUser)
                self.centcmod.update(util.centcmod.CORRPATTERN,self.corrPattern)

    def takeCentWeight(self,control):
        """If centWeight is a string (eg make or something similar), use a default SH spot pattern as the reference centroid weighting.
        """
        if type(self.centWeight)==type(""):
            #Need to create the centroid weighting.
            if self.corrPattern==None:
                c=control["cal_source"]
                control["cal_source"]=1
                steps=self.linearSteps
                self.linearSteps=None
                if self.centcmod==None:
                    raise Exception("Not yet sorted for non-cmod type things")
                self.reorderedPhs[:]=0
                self.runCalc(control)
                control["cal_source"]=c
                self.linearSteps=steps
                self.centWeight=self.cmodbimg.copy()
                self.centcmod.update(util.centcmod.CENTWEIGHT,self.centWeight)

def compute(phase,nsubx):
    """Simple interface..."""
    c=centroid(nsubx)
    c.run(phase)
    return c

def createAiryDisc(npup,halfwidth=1.,xoff=0.,yoff=0.):
    """Creates an airy disk pattern central on an array npup pixels wide,
    with a width (2*halfwidth) specified in pixels.
    xoff/yoff can be used to shift the airy disc around the array:
    central by default for both odd and even sized arrays.
    This can be used for convolving with the high light level SHS spots
    (ie before noise addition).
    """
    import scipy.special
    #first zero of j1 occurs at 3.8317 (http://mathworld.wolfram.com/BesselFunctionZeros.html).
    #However, this should occur at 1.22lambda/D.  So need to scale this.
    center=(npup%2-1.)/2
    dist=util.dist.dist(npup,dy=center+yoff,dx=center+xoff)*3.8317/halfwidth
    disc=scipy.special.j1(dist)/dist
    disc=numpy.where(dist==0,0.5,disc)**2#remove NaN and square...
    return disc


class fpgaCentStateInformation:
    def __init__(self):
        self.calSource=None
        self.seedLoaded=None
        self.pupil=None
        self.nimg=None
        self.fftsize=None
        self.cenmask=None
        self.doPartialFPGA=None
        self.fittedSubaps=None
        self.nIntegrations=None
        self.nsubx=None
        self.phasesize=None
        self.sigFPGA=None
        self.fpgaSkybrightness=None
        self.noiseFloor=None
        self.fpga_readbg=None
        self.fpga_readnoise=None
        self.symtype=None

class CalData:
    def __init__(self,xc,xr,yc,yr):
        self.xc=xc#x centroids
        self.xr=xr#x range
        self.yc=yc
        self.yr=yr
        self.indx=[]
        self.xindx=None
        self.yindx=None
