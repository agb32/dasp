import cmod.fft
import util.arrayFromArray
import numpy
from util.flip import fliparray2
from util.dist import dist
import time
import util.FITS
import os
import tempfile
import cmod.binimg
haveFPGA=1
try:
    import fpga
except:
    print "FPGA module not installed."
    haveFPGA=0

class science:
    """A class of utilities related to science imaging.  This can be used
    stand-alone or as part of the aosim/science/science.py module.
    realPup - used if phase is binned before image creation.
    Binning phase can be used to reduce computation time, or eg make it fit into the FPGA... for a 42m telescope, 1024x1024 phase values are fine, binning down from 2048 or 4096 gives virtually no difference to strehl ratio.  This is because at 1024, each phase size is about 4cm, which is less than half r0... at 2048, this is even smaller, so doesn't change anything.
    InboxDiamList is a list of diameters (arcsec) for which the ensquared energy should be calculated.
    phaseMultiplier is used for cases where phase wavelength is different from the wavelength that this science object is at.

    """
    def __init__(self,npup,nfft,pup,nimg=None,tstep=0.005,atmosPhaseType="phaseonly",keepDiffPsf=0,pix_scale=1.,fitsFilename=None,diffPsfFilename=None,scinSamp=1,sciPSFSamp=1,scienceListsSize=128,debug=None,timing=0,allocateMem=1,realPup=None,useFPGA=0,waitFPGA=0,waitFPGATime=0,fpDataType="f",calcRMSFile=None,inboxDiamList=[0.2],sciFilename=None,saveFileString=None,nthreads=1,histFilename=None,phaseMultiplier=1,luckyNSampFrames=1,luckyFilename=None,luckyImgFilename=None,luckyImgSize=None,luckyHistorySize=10):
        self.fftPlan=None
        self.nfft=nfft
        self.npup=npup
        if nimg==None:
            nimg=nfft
        self.nimg=nimg
        if (self.nfft%self.nimg)!=0:
            raise Exception("WARNING: util.sci.py nfft not a multiple of nimg")
        self.phaseTilt=None
        self.fpgaUserPrecision=4
        if npup>=512:
            self.fpgaUserPrecision-=1
        if npup>=1024:
            self.fpgaUserPrecision-=1
        self.pup=pup
        self.tstep=tstep
        self.atmosPhaseType=atmosPhaseType
        if self.atmosPhaseType not in ["phaseonly","phaseamp","realimag"]:
            raise Exception("science: atmosPhaseType not known %s"%self.atmosPhaseType)
        self.keepDiffPsf=keepDiffPsf
        self.pix_scale=pix_scale
        self.fitsFilename=fitsFilename#these are used when the science.science object is destroyed.
        self.sciFilename=sciFilename
        self.diffPsfFilename=diffPsfFilename
        self.scinSamp=scinSamp
        self.sciPSFSamp=sciPSFSamp
        self.computeOTF=0#compute OTF for strehl calcs
        self.historyListsSize=scienceListsSize
        self.luckyHistorySize=luckyHistorySize
        self.inboxDiamList=inboxDiamList
        self.saveFileString=saveFileString
        if type(self.saveFileString)!=type(""):
            self.saveFileString=""
        self.debug=debug
        self.nthreads=nthreads
        self.timing=timing
        if type(realPup)!=type(None):#this is used if the phase is binned before image creation.
            self.realPupBinned=numpy.zeros((self.npup,self.npup),numpy.float32)
            cmod.binimg.binimg(realPup.astype("f"),self.realPupBinned)
            self.realPupBinned[:,]=numpy.where(self.realPupBinned==0,1,self.realPupBinned)
        self.useFPGA=useFPGA
        self.waitFPGA=waitFPGA
        self.waitFPGATime=waitFPGATime
        self.fpgaWaitCycles=0#number of cycles waited for FPGA for - check occasionalyl to make sure we;re not waiting too long...
        self.fpDataType=fpDataType#floating point type
        self.cpDataType=fpDataType.upper()#complex type
        if (self.nfft/self.nimg)%2==0 and self.nfft!=self.nimg:#even binning
            self.computePhaseTilt()

        self.calcRMSFile=calcRMSFile
        self.phaseRMS=0
        self.phaseRMSsum=0.
        self.phaseRMSsum2=0.
        self.phaseRMScnt=0
        self.pupsum=numpy.sum(numpy.sum(self.pup.astype("l")))        
        self.idxPup=numpy.nonzero(self.pup.flat) #Flat version of the pixels belonging to the pupil#
        self.initRadialProfile()
        self.initEnsquaredEnergyProfile()
        #this.inbox=Numeric.zeros(nfft/2,Numeric.Float64)# In box energy
        self.wid=self.rSquare*self.pix_scale #Square aperture size in arcsec. X Axis of inbox array, used for plotting		
        #self.strehl_hist=[] ## History of Strehl ratio (for plotting)
        #self.d50_hist=[] ##History of d50
        #self.fwhm_hist=[] ##History of fwhm
        #self.inbox_hist=[] ##History of inbox_2
        self.histFilename=histFilename
        self.history=None#history of the science parameters...
        self.dictScience={}
        self.n_integn=0
        self.isamp=0
        self.psfSamp=0
        self.PSFTime=0.
        self.phaseMultiplier=phaseMultiplier
        self.longExpPSF=numpy.zeros((nimg,nimg),self.fpDataType)#was float64
        self.longExpImg=numpy.zeros((nimg,nimg),self.fpDataType)# Long exposure image array (cannot be shared)
        self.luckyImg=None
        self.luckyCnt=0
        self.luckyNSampFrames=luckyNSampFrames#Note, every sciPSFSamp (which defaults to 1, not scinsamp) frames, the lucky image will be integrated, until this has been done luckyNSampFrames, at which point results are calculated and the luckyImg zeroed.
        self.luckyDict={}
        self.luckyHistory=None
        self.luckyFilename=luckyFilename
        self.luckyImgFilename=luckyImgFilename
        self.luckyImgSize=luckyImgSize
        self.luckyFile=None
        if self.useFPGA:
            self.testFPGAUsage()
        else:
            self.canUseFPGA=0
        if allocateMem:
            self.initMem()
            self.initProfiles()
        self.longExpPSFView=None
        self.instImgView=None

    def __del__(self):
        if self.fftPlan!=None:
            cmod.fft.FreePlan(self.fftPlan)
            cmod.fft.CleanUp()

    def computePhaseTilt(self):
        """A phase tilt is necessary since binning with an even number...
        This ensures the PSF is squarely spaced before binning"""
        bf=self.nfft/self.nimg
        tmp=float(self.npup)/float(self.nfft)*numpy.pi*(bf-1)
        xtiltfn=((numpy.fromfunction(lambda x,y:y,(self.npup,self.npup))-float(self.npup)/2.+0.5)/float(self.npup)).astype(self.fpDataType)# subap tilt fn
        self.phaseTilt=(self.pup*tmp*(xtiltfn+numpy.transpose(xtiltfn))).astype(self.fpDataType)

        
    def initMem(self,fftTmp=None,pupilAmplitude=None,focusAmplitude=None,tempImg=None,instImg=None,phs=None,binImg=None):
        """Initialise memories needed... here, if memories are passed in,
        these are used (ie when sharing resources), otherwise, new memories
        are allocated."""
        nfft=self.nfft
        npup=self.npup
        nimg=self.nimg
        if self.canUseFPGA:#pupilAmplitude should be the fpga accessible memory.
            if self.atmosPhaseType=="phaseonly":
                atmosfactor=1#phaseonly...
            else:
                raise Exception("util.science: atmosphasetype: not phaseonly")
            
        if type(fftTmp)==type(None):#scratch area for acml FFT routine...
            self.fftTmp=numpy.zeros((nfft**2+10*nfft,),self.cpDataType)#Numeric.Complex64
        elif fftTmp.shape!=(nfft**2+10*nfft,):
            self.fftTmp=util.arrayFromArray.arrayFromArray(fftTmp,(nfft**2+10*nfft,),self.cpDataType)#Numeric.Complex64
        else:
            self.fftTmp=fftTmp
        if type(pupilAmplitude)==type(None):
            self.pupilAmplitude=numpy.zeros((nfft,nfft),self.cpDataType)#Numeric.Complex64
        elif pupilAmplitude.shape!=(nfft,nfft):
            self.pupilAmplitude=util.arrayFromArray.arrayFromArray(pupilAmplitude,(nfft,nfft),self.cpDataType)#Numeric.Complex64
        else:
            self.pupilAmplitude=pupilAmplitude
        #self.pupilAmplitude.savespace(1)
        if focusAmplitude==None:#an in place FFT will be done (slower)
            self.focusAmplitude=self.pupilAmplitude
        elif focusAmplitude.shape!=(nfft,nfft):
            self.focusAmplitude=util.arrayFromArray.arrayFromArray(focusAmplitude,(nfft,nfft),self.cpDataType)#Numeric.Complex64
        else:
            self.focusAmplitude=focusAmplitude

            

        #initialise the FFT (fftw module)
        cmod.fft.InitialiseThreading(self.nthreads)
        self.fftPlan=cmod.fft.Plan(self.pupilAmplitude,self.focusAmplitude)
        
        if type(tempImg)==type(None):
            self.tempImg=numpy.zeros((nfft,nfft),self.fpDataType)#was float64
        elif tempImg.shape!=(nfft,nfft):
            self.tempImg=util.arrayFromArray.arrayFromArray(tempImg,(nfft,nfft),self.fpDataType)#was float64
        else:
            self.tempImg=tempImg
        if self.nfft==self.nimg:
            self.binImg=self.tempImg
        else:
            if type(binImg)==type(None):
                self.binImg=numpy.zeros((self.nimg,self.nimg),self.fpDataType)
            elif binImg.shape!=(self.nimg,self.nimg):
                self.binImg=util.arrayFromArray.arrayFromArray(binImg,(self.nimg,self.nimg),self.fpDataType)
            else:
                self.binImg=binImg
        
        instImg=numpy.array(instImg)#xxx
        if type(instImg)==type(None):
            self.instImg=numpy.zeros((self.nimg,self.nimg),self.fpDataType)#was Numeric.Float64
        else:
            self.instImg=util.arrayFromArray.arrayFromArray(instImg,(nimg,nimg),self.fpDataType)#Numeric.Float64
    
        if type(phs)==type(None):
            if self.atmosPhaseType=="phaseonly":
                self.phs=numpy.zeros((npup,npup),self.fpDataType)#Numeric.Float64
            elif self.atmosPhaseType=="phaseamp":
                self.phs=numpy.zeros((2,npup,npup),self.fpDataType)#Numeric.Float64
            elif self.atmosPhaseType=="realimag":
                self.phs=numpy.zeros((npup,npup),self.cpDataType)#Numeric.Complex64
        else:
            shape={"phaseonly":(npup,npup),"phaseamp":(2,npup,npup),"realimag":(npup,npup)}[self.atmosPhaseType]
            dtype={"phaseonly":self.fpDataType,"phaseamp":self.fpDataType,"realimag":self.cpDataType}[self.atmosPhaseType]
            phs=numpy.array(phs)#xxx
            self.phs=util.arrayFromArray.arrayFromArray(phs,shape,dtype)
        self.phsInt=util.arrayFromArray.arrayFromArray(self.phs,(npup,npup),numpy.int32)
##         if type(longExpPSF)==type(None):
##             self.longExpPSF=Numeric.zeros((nfft,nfft),Numeric.Float64)
##         elif longExpPSF.shape!=(nfft,nfft):
##             self.longExpPSF=cmod.utils.arrayFromArray(longExpPSF,(nfft,nfft),Numeric.Float64)
##         else:
##             self.longExpPSF=longExpPSF
        if self.canUseFPGA:
            self.fpgaPupilMask=(0x7f800000*(1-self.pup)).astype(numpy.int32)
    def initialiseFPGA(self,fpid=None,fpgaInfo=None,ignoreFailure=0,fpgaBitFile=None):
        self.fpid=fpid
        self.fpgaInfo=fpgaInfo
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
                        print "Warning: util.science: failed to initialise FPGA"
                    else:
                        print "Failed to open FPGA device..."
                        print os.environ
                        raise
                if self.fpgaBinaryLoaded:
                    fpga.load(self.fpid,fpgaBitFile)
                    fpga.reset(self.fpid)
                    time.sleep(0.001)
                    fpga.start(self.fpid)
                    time.sleep(0.001)
                    fpga.writeReg(self.fpid,0x2,6)#stop the fpga pipeline
                    self.fpgaInfo=fpgaSciStateInformation()
            else:
                self.fpgaBinaryLoaded=1
        if self.fpgaBinaryLoaded==0:
            self.fpid=None
            self.fpgaInfo=None
            self.canUseFPGA=0
        return self.fpid,self.fpgaInfo
    def testFPGAUsage(self):
        self.canUseFPGA=haveFPGA
        if self.nfft not in [8,16,32,64,128,256,512,1024,2048]:
            print "WARNING: cannot use FPGA - illegal nfft value"
            self.canUseFPGA=0
        if self.nimg!=self.nfft:
            print "WARNING: cannot use FPGA - illegal nimg value"
            self.canUseFPGA=0
        if self.npup>1024:
            print "WARNING: cannot use FPGA - npup>1024"
            self.canUseFPGA=0
        if self.fpDataType!=numpy.float32:
            print "WARNING: cannot use FPGA - must be float32"
            self.canUseFPGA=0
        #if self.cpDataType!=numpy.Complex32:
        #    self.canUseFPGA=0
        if self.atmosPhaseType!="phaseonly":
            print "WARNING: cannot use FPGA - atmosPhaseType must be phaseonly"
            self.canUseFPGA=0
        return self.canUseFPGA
    def setupFPGAArray(self,fpid,fftsize):
        """Set up FPGA accessible region.  The phase read into the fpga can be overwritten by the image..."""
        if fftsize>2048:
            raise Exception("fft size too large for science fpga (2048 max, requested %d)"%fftsize)
        fpgaarr=fpga.mallocHostMem(fpid,(fftsize**2*4,),numpy.int8)
        return fpgaarr
    def fpgaLoadRegs(self):
        fpga.writeReg(self.fpid,2,6)#stop the pipe
        fpga.writeReg(self.fpid,self.npup,512)
        fpga.writeReg(self.fpid,self.nfft,513)
        fpga.writeReg(self.fpid,self.fpgaUserPrecision,514)
        self.fpgaInfo.npup=self.npup
        self.fpgaInfo.nfft=self.nfft
        self.fpgaInfo.userPrecision=self.fpgaUserPrecision
        fpga.writeAddr(self.fpid,self.phs,1)#input address
        fpga.writeReg(self.fpid,self.npup*self.npup/2,2)#read into fpga size
        fpga.writeAddr(self.fpid,self.instImg,3)#output addr
        fpga.writeReg(self.fpid,self.nfft*self.nfft/2,4)#output size
        fpga.writeReg(self.fpid,64,6)#reset the fifos.
        
    def initProfiles(self,useFPGA=0):
        """This should be called after the memory has been set up..."""
        useFPGA=useFPGA and self.canUseFPGA
        self.diffPSF=self.computeDiffPSF(useFPGA)
        self.diffPsfRadialProfile=self.computeRadialProfileAndEncircledEnergy(self.diffPSF)[0,:,]
        if self.diffPsfFilename!=None:#save the diffraction limited PSF.
            util.FITS.Write(self.diffPSF,self.diffPsfFilename)
        #self.dlPsfRadialProfile=self.computeRadialProfileAndEncircledEnergy(self.diffPSF)[0,:,]
        self.diffn_core_en=float(self.diffPSF[self.nimg/2,self.nimg/2])
        self.diffOTFSum=numpy.fft.fft2(numpy.fft.fftshift(self.diffPSF),s=(self.diffPSF.shape[0]*2,self.diffPSF.shape[1]*2)).sum()
        print "Diffraction OTF sum: %s"%str(self.diffOTFSum)
        print "todo - return from computeEnsquaredEnergy - allocate once"
        self.diffPsfEnsquaredProfile=self.computeEnsquaredEnergy(self.diffPSF)
        #perform an acml FFT initialisation.
        #Actually not needed for inplace_fft2d...
##     def initProfiles(self,useFPGA=0):
##         """This should be called after the memory has been set up..."""
##         useFPGA=useFPGA and self.canUseFPGA
##         self.diffPSF=self.computeDiffPSF(useFPGA)
##         self.diffPsfRadialProfile=self.computeRadialProfileAndEncircledEnergy(self.diffPSF)[0,:,]
##         if self.diffPsfFilename!=None:#save the diffraction limited PSF.
##             util.FITS.Write(self.diffPSF,self.diffPsfFilename)
##         #self.dlPsfRadialProfile=self.computeRadialProfileAndEncircledEnergy(self.diffPSF)[0,:,]
##         self.diffn_core_en=float(self.diffPSF[self.nfft/2,self.nfft/2])
##         print "todo - return from computeEnsquaredEnergy - allocate once"
##         self.diffPsfEnsquaredProfile=self.computeEnsquaredEnergy(self.diffPSF)
##         #perform an acml FFT initialisation.
##         #Actually not needed for inplace_fft2d...
        

    def initRadialProfile(self):
        """Creates and initialises arrays to compute the radial
        profile and the encircled energy profile of the PSF Must be
        first called before calling
        computeRadialProfileAndEncircledEnergy
        """	
        ####  We first create a map of the square of the distances to the center of the PSF image
        ####  The map of distances is computed with the dist function (see aosim/utils/dist.py for help)
        ####  We use the square because we have then integer indexes
        r2=dist(self.nimg,natype=self.fpDataType,sqrt=0).ravel()
        #r2*=r2 ##no longer sqrts in dist... square : we have integer numbers
        self.nnRad=numpy.argsort(r2) ##we flatten the grid of distances and sort the values
        r2r=numpy.take(r2,self.nnRad) ## we sort the grid of distances
        self.xRad=numpy.nonzero(difYorick(r2r))[0] ##we look when the grid of distances change of value
        #print self.xRad,type(self.xRad),type(self.xRad[0])
        self.tabNbPointsRad=difYorick(self.xRad) ## number of points per bin
        self.rRad=numpy.take(numpy.sqrt(r2r),self.xRad) ##radius (in pixels) giving the horizontal axis for radial and encircled energy profiles
        self.rad=self.rRad*self.pix_scale

##     def initRadialProfile(self):
##         """Creates and initialises arrays to compute the radial
##         profile and the encircled energy profile of the PSF Must be
##         first called before calling
##         computeRadialProfileAndEncircledEnergy
##         """	
##         ####  We first create a map of the square of the distances to the center of the PSF image
##         ####  The map of distances is computed with the dist function (see aosim/utils/dist.py for help)
##         ####  We use the square because we have then integer indexes
##         r2=dist(self.nfft,natype=self.fpDataType,sqrt=0)
##         r2=Numeric.array(r2)
##         #r2*=r2 ##no longer sqrts in dist... square : we have integer numbers
##         self.nnRad=Numeric.argsort(r2.flat) ##we flatten the grid of distances and sort the values
##         r2r=Numeric.take(r2.flat,self.nnRad) ## we sort the grid of distances
##         self.xRad=Numeric.nonzero(dif_(r2r)) ##we look when the grid of distances change of value
##         self.tabNbPointsRad=dif_(self.xRad) ## number of points per bin
##         self.rRad=Numeric.take(Numeric.sqrt(r2r),self.xRad) ##radius (in pixels) giving the horizontal axis for radial and encircled energy profiles
##         self.rad=self.rRad*self.pix_scale

    def computeRadialProfileAndEncircledEnergy(self,psf):
        """Computes radial end encircled energy profiles
            Inputs :
             - psf : the input whose profiles are wanted

            Outputs :
              - result : the array storing profiles
                the first line, ie result[0,] stores the radial profile
                the second line, ie result[1,] stores the encircled energy profile
        """
        #### We flatten the PSF and sort the values
        #psf=numpy.array(psf)
        #print type(psf)
        psfSort=numpy.take(psf.ravel(),self.nnRad)
        #### We compute the encircled energy of the PSF
        tabEncircledEnergy=numpy.take(numpy.cumsum(psfSort),self.xRad)
        #### We compute the radial profile
        tabEnergyPerPixBin=difYorick(tabEncircledEnergy) ##to have the energy per pixel bin
        #tabEnergyPerPixBin.savespace(1)#prevent conversion to double.
        profil=tabEnergyPerPixBin/self.tabNbPointsRad

        #### Allocation of the return result
        result=numpy.zeros((2,len(tabEncircledEnergy)),dtype=psf.dtype)#savespace=1)

        #### We first store the radial profile in line 1 of the returned array
        result[0,0]=psfSort[0]
        result[0,1:,]=profil[:,] ##radial profile of PSF

        #### We  store the encircled energy profile in line 2 of the returned array
        result[1,:,]=tabEncircledEnergy[:,]
        return result

##     def computeRadialProfileAndEncircledEnergy(self,psf):
##         """Computes radial end encircled energy profiles
##             Inputs :
##              - psf : the input whose profiles are wanted

##             Outputs :
##               - result : the array storing profiles
##                 the first line, ie result[0,] stores the radial profile
##                 the second line, ie result[1,] stores the encircled energy profile
##         """
##         #### We flatten the PSF and sort the values
##         psfSort=Numeric.take(psf.flat,self.nnRad)
##         #### We compute the encircled energy of the PSF
##         tabEncircledEnergy=Numeric.take(Numeric.cumsum(psfSort),self.xRad)
##         #### We compute the radial profile
##         tabEnergyPerPixBin=dif_(tabEncircledEnergy) ##to have the energy per pixel bin
##         tabEnergyPerPixBin.savespace(1)#prevent conversion to double.
##         profil=tabEnergyPerPixBin/self.tabNbPointsRad

##         #### Allocation of the return result
##         result=Numeric.zeros((2,len(tabEncircledEnergy)),typecode=psf.typecode(),savespace=1)
##         print type(profil),profil.shape,result.shape,profil.dtype.char,result.typecode()
##         #### We first store the radial profile in line 1 of the returned array
##         result[0,0]=psfSort[0]
##         result[0,1:,]=profil[:,].astype("f") ##radial profile of PSF

##         #### We  store the encircled energy profile in line 2 of the returned array
##         result[1,:,]=tabEncircledEnergy[:,]
##         return result

    def initEnsquaredEnergyProfile(self):
        """Creates and initialises arrays to fastly compute the ensquared energy profile
        Must be called before calling computeEnsquaredEnergy
        """	
        ####  We first create a pixel map of concentric square apertures
        tabx=numpy.arange(self.nimg)-self.nimg/2;
        r2=numpy.maximum(numpy.absolute(tabx[:,numpy.newaxis]),numpy.absolute(tabx[numpy.newaxis,:,]))

        ##we flatten the grid of distances and sort the values
        self.nnSquare=numpy.argsort(r2.ravel())
        r2r=numpy.take(r2.ravel(),self.nnSquare) ## we sort the grid of distances
        self.xSquare=numpy.nonzero(difYorick(r2r))[0] ##we look when the grid of distances change of value
        self.rSquare=numpy.take(r2r,self.xSquare)*2+1 ##aperture size (in pixels) giving the horizontal axis for the ensquared energy profile
    
##     def initEnsquaredEnergyProfile(self):
##         """Creates and initialises arrays to fastly compute the ensquared energy profile
##         Must be called before calling computeEnsquaredEnergy
##         """	
##         ####  We first create a pixel map of concentric square apertures
##         tabx=Numeric.arange(self.nfft)-self.nfft/2;
##         r2=Numeric.maximum(Numeric.absolute(tabx[:,Numeric.NewAxis]),Numeric.absolute(tabx[Numeric.NewAxis,:,]))

##         ##we flatten the grid of distances and sort the values
##         self.nnSquare=Numeric.argsort(r2.flat)
##         r2r=Numeric.take(r2.flat,self.nnSquare) ## we sort the grid of distances
##         self.xSquare=Numeric.nonzero(dif_(r2r)) ##we look when the grid of distances change of value
##         self.rSquare=Numeric.take(r2r,self.xSquare)*2+1 ##aperture size (in pixels) giving the horizontal axis for the ensquared energy profile
    

    def computeEnsquaredEnergy(self,psf):
        """Computes ensquared energy profile
        """
        #### We flatten the diffraction limited PSF and sort the values
        psfSort=numpy.take(psf.ravel(),self.nnSquare)
        #### We compute the encircled energy of the diffraction limited PSF 
        tabEnsquaredEnergy=numpy.take(numpy.cumsum(psfSort),self.xSquare)
        #### We convert the profile into psf type and return the array
        result=tabEnsquaredEnergy.astype(psf.dtype)
        return result
##     def computeEnsquaredEnergy(self,psf):
##         """Computes ensquared energy profile
##         """
##         #### We flatten the diffraction limited PSF and sort the values
##         #psf=Numeric.array(psf)
##         psfSort=Numeric.take(psf.flat,self.nnSquare)

##         #### We compute the encircled energy of the diffraction limited PSF 
##         tabEnsquaredEnergy=Numeric.take(Numeric.cumsum(psfSort),self.xSquare)

##         #### We convert the profile into psf type and return the array
##         result=tabEnsquaredEnergy.astype(psf.typecode())
##         return result
    def calcRMS(self,phase,mask=None):
        """calcs rms in radians (or whatever phase is in)"""
        if type(mask)==type(None):
            rms=numpy.sqrt(numpy.average(phase.flat**2)-numpy.average(phase.flat)**2)
        else:
            p=(phase*mask).ravel()
            p2=numpy.sum(p*p)
            s=numpy.sum(mask)
            m=numpy.sum(p)/s
            rms=numpy.sqrt(p2/s-m*m)
        return rms

##     def calcRMS(self,phase,mask=None):
##         """calcs rms in radians (or whatever phase is in)"""
##         if type(mask)==type(None):
##             rms=Numeric.sqrt(Numeric.average(phase.flat**2)-Numeric.average(phase.flat)**2)
##         else:
##             p=(phase*mask).flat
##             p2=Numeric.sum(p*p)
##             s=Numeric.sum(mask.astype("d").flat)
##             m=Numeric.sum(p)/s
##             rms=Numeric.sqrt(p2/s-m*m)
##         return rms

    def computeDiffPSF(self,useFPGA=0):
        """computeDiffPSF function : computes the diffraction limited PSF
        Function written  by FA, heavily modified by AGB"""
        useFPGA=(useFPGA and self.canUseFPGA)
        if self.atmosPhaseType=="phaseamp":
            self.computeShortExposurePSF(numpy.array([1,1],self.fpDataType),useFPGA)#places result into instImg
        else:
            self.computeShortExposurePSF(numpy.array([1],self.fpDataType),useFPGA)
        if self.keepDiffPsf:
            self.diffPSF=self.instImg.copy()
        else:
            self.diffPSF=self.instImg#will be overwritten next time computeShortExposurePSF is called...
        return self.diffPSF
##     def computeDiffPSF(self,useFPGA=0):
##         """computeDiffPSF function : computes the diffraction limited PSF
##         Function written  by FA, heavily modified by AGB"""
##         useFPGA=useFPGA and self.canUseFPGA
##         if self.atmosPhaseType=="phaseamp":
##             self.computeShortExposurePSF(Numeric.array([1,1],self.fpDataType),useFPGA)#places result into instImg
##         else:
##             self.computeShortExposurePSF(Numeric.array([1],self.fpDataType),useFPGA)
##         if self.keepDiffPsf:
##             self.diffPSF=self.instImg.copy()
##         else:
##             self.diffPSF=self.instImg#will be overwritten next time computeShortExposurePSF is called...
##         return self.diffPSF


    def computeShortExposurePSF(self,phs,useFPGA=0):
        """computeShortExposurePSF function : computes short exposure AO corrected PSF
        Modifications made by FA"""

##        ##We copy the arrays to have contigous arrays and convert it into float64
##        phs2=phs.copy()
##        phs2=phs2.astype("d")
##        self.phs2=phs2
##        pup2=self.pup.copy()
##        pup2=pup2.astype("d")
##        self.pup2=pup2
        t1=time.time()
        npup=self.npup ##to have the dimensions of the input phase array
        if self.phaseTilt!=None:#a phase tilt is needed if binning by an even factor (to put the central spot in a single pixel)
            phs=phs+self.phaseTilt
        
        # print self.pupilAmplitude.shape,self.phs.shape,self.pup.shape
        # We fill the complex amplitude
        if self.atmosPhaseType=="phaseonly":
            self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
            self.pupilAmplitude[:npup,npup:,]=0.
            self.pupilAmplitude.real[:npup,:npup]=self.pup*numpy.cos(phs)
            self.pupilAmplitude.imag[:npup,:npup]=self.pup*numpy.sin(phs)
        elif self.atmosPhaseType=="phaseamp":#phs[1] is amplitude, phs[0] is phase
            self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
            self.pupilAmplitude[:npup,npup:,]=0.
            self.pupilAmplitude.real[:npup,:npup]=self.pup*numpy.cos(phs[0])*phs[1]
            self.pupilAmplitude.imag[:npup,:npup]=self.pup*numpy.sin(phs[0])*phs[1]
        elif self.atmosPhaseType=="realimag":#phs in real/imag already.
            self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
            self.pupilAmplitude[:npup,npup:,]=0.
            self.pupilAmplitude[:npup,:npup]=self.pup*phs


        ##We call the FFTW function
        cmod.fft.ExecutePlan(self.fftPlan)#fft from pupilAmplitude to focusAmplitude
        #util.acmlfft.inplace_fft2d(self.pupilAmplitude,self.fftTmp)
        #nrfft.fftw_execute(self.fftPlan)
        ##We compute the intensity in the focal plane
        self.tempImg[:,]=(numpy.absolute(self.focusAmplitude))#**2
        self.tempImg[:,]*=self.tempImg[:,]#now square it (slightly faster)
        ##print phs2.shape,phs2.typecode(),self.tempImg.shape,self.tempImg.typecode(),pup2.shape,pup2.typecode()
        ##We compute the PSF by using fliparray2
        #self.instImg=fliparray2(self.tempImg)# Flip quadrants with definition consistent with FFT coordinate definition
        if self.nimg!=self.nfft:#bin the image...
            cmod.binimg.binimg(self.tempImg,self.binImg)
        fliparray2(self.binImg,self.instImg)# Flip quadrants with definition consistent with FFT coordinate definition
        self.instImg/=numpy.sum(self.instImg) ##We normalise the instantaneous PSF to 1
        t2=time.time()
        self.PSFTime=t2-t1
        #print numpy.sum(numpy.sum(self.instImg))
        return self.instImg
                
##     def computeShortExposurePSF(self,phs,useFPGA=0):
##         """computeShortExposurePSF function : computes short exposure AO corrected PSF
##         Modifications made by FA"""

## ##        ##We copy the arrays to have contigous arrays and convert it into Float64
## ##        phs2=phs.copy()
## ##        phs2=phs2.astype("d")
## ##        self.phs2=phs2
## ##        pup2=self.pup.copy()
## ##        pup2=pup2.astype("d")
## ##        self.pup2=pup2
##         t1=time.time()
##         npup=self.npup ##to have the dimensions of the input phase array
##         if useFPGA and self.useFPGA:#mask the unused parts of the pupil.
##             if phs is not self.phs:# check that the phase is in fpga array
##                 print "copying phase array"
##                 self.phs[:,]=phs
##             # now mask out the pupil 
##             self.phsInt|=self.fpgaPupilMask#phsInt points to same memory as phs
##             fpga.writeReg(self.fpid,2,6)#stop pipe
##             fpga.writeReg(self.fpid,64,6)#reset fifos
##             if self.npup!=self.fpgaInfo.npup or self.nfft!=self.fpgaInfo.nfft or self.fpgaUserPrecision!=self.fpgaInfo.userPrecision:
##                 self.fpgaLoadRegs()
##             fpga.writeReg(self.fpid,1,6)#set it going.
##             if self.waitFPGA:
##                 self.fpgaWaitCycles=0
##                 time.sleep(self.waitFPGATime)
##                 v=fpga.readReg(self.fpid,5)
##                 while v!=7:
##                     v=fpga.readReg(self.fpid,5)
##                     #print hex(v)
##                     self.fpgaWaitCycles+=1
##                 if fpga.readReg(self.fpid,515)==1:#overflow error... try again.
##                     print "FPGA overflow error, reducing precision (science)"
##                     self.fpgaUserPrecision-=1
##                     if self.fpgaUserPrecision<0:
##                         self.fpgaUserPrecision=0
##                         raise Exception("FPGA science FFT precision overflow")
##                     fpga.reset(self.fpid)#clear the overflow bit.
##                     fpga.start(self.fpid)
##                     self.fpgaLoadRegs()
##                     self.computeShortExposurePSF(self.phs,useFPGA)#re run the calc.
##         else:
##             # print self.pupilAmplitude.shape,self.phs.shape,self.pup.shape
##             # We fill the complex amplitude
##             if self.atmosPhaseType=="phaseonly":
##                 self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
##                 self.pupilAmplitude[:npup,npup:,]=0.
##                 self.pupilAmplitude.real[:npup,:npup]=self.pup*Numeric.cos(phs)
##                 self.pupilAmplitude.imag[:npup,:npup]=self.pup*Numeric.sin(phs)
##             elif self.atmosPhaseType=="phaseamp":#phs[1] is amplitude, phs[0] is phase
##                 self.pupilAmplitude[
## npup:,]=0.#clear array (may get changed by fft)
##                 self.pupilAmplitude[:npup,npup:,]=0.
##                 self.pupilAmplitude.real[:npup,:npup]=self.pup*Numeric.cos(phs[0])*phs[1]
##                 self.pupilAmplitude.imag[:npup,:npup]=self.pup*Numeric.sin(phs[0])*phs[1]
##             elif self.atmosPhaseType=="realimag":#phs in real/imag already.
##                 self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
##                 self.pupilAmplitude[:npup,npup:,]=0.
##                 self.pupilAmplitude[:npup,:npup]=self.pup*phs


##             ##We call the FFTW function
##             util.acmlfft.inplace_fft2d(self.pupilAmplitude,self.fftTmp)
##             #nrfft.fftw_execute(self.fftPlan)
##             ##We compute the intensity in the focal plane
##             self.tempImg[:,]=(Numeric.absolute(self.focusAmplitude))**2
##             ##print phs2.shape,phs2.typecode(),self.tempImg.shape,self.tempImg.typecode(),pup2.shape,pup2.typecode()
##             ##We compute the PSF by using fliparray2
##             #self.instImg=fliparray2(self.tempImg)# Flip quadrants with definition consistent with FFT coordinate definition
##             fliparray2(self.tempImg,self.instImg)# Flip quadrants with definition consistent with FFT coordinate definition
##             self.instImg/=Numeric.sum(Numeric.sum(self.instImg)) ##We normalise the instantaneous PSF to 1
##         t2=time.time()
##         self.PSFTime=t2-t1
##         #print Numeric.sum(Numeric.sum(self.instImg))
##         return self.instImg


    def computeScientificParameters(self,longExpPSF=None,dictScience=None):
        """computeScientificParameters function : computes FWHM, 
           Modifications made by FA
           The longExpPSF should be normalised to 1, i.e.sum()==1.
           """
        if longExpPSF==None:
            longExpPSF=self.longExpPSF
        if dictScience==None:
            dictScience=self.dictScience
        #window(2,wait=1)
        #fma()
        #pli(self.longExpPSF)
                
        ## We calculate PSF parameters ###
        ## We start by the Strehl Ratio
        nfft=self.nfft
        nimg=self.nimg
        strehl=longExpPSF[nimg/2,nimg/2]/self.diffn_core_en
        dictScience['strehl']=strehl
        ##print "Strehl=%g"%(strehl)
        dictScience['strehlPeak']=max(longExpPSF.flat)/self.diffn_core_en
        if self.computeOTF:#compute Strehl by OTF method.
            sl=numpy.fft.fftshift(longExpPSF)
            dictScience['strehlOTF']=numpy.fft.fft2(sl,s=(sl.shape[0]*2,sl.shape[1]*2)).sum()/self.diffOTFSum

        ## We now compute the FWHM
        ## First we compute the radial profile
        ##Pbs with interp function : we have to compute it manually :-((
        profiles=self.computeRadialProfileAndEncircledEnergy(longExpPSF)
        x=self.rad[:,]
        y=profiles[0,:,]
        ##we look for the index of x just before and after y[0]/2
        t=numpy.greater(y,y[0]/2)
        try:
            i1=numpy.nonzero(t)[0][-1]
        except:
            print "ERROR in util/science.py - no non-zero elements... Oh well!"
            i1=0
        x1=x[i1]
        t=numpy.less(y,y[0]/2)
        try:
            i2=numpy.nonzero(t)[0][0]
        except:
            print "ERROR2 in util/science.py - no non-zero elements... Oh well! "
            i2=0
        x2=x[i2]
        ##we compute the equation of the line between those two points
        if x2==x1:
            x2+=1
        #print y[i2].shape,y[i1].shape,type(y),y.shape,i2,i1
        a=(y[i2]-y[i1])/(x2-x1)
        b=y[i2]-a*x2
        ##we then compute the FWHM
        fwhm=2*(y[0]/2-b)/a
        ##print "FWHM=%g"%(fwhm)
        dictScience['FWHM']=fwhm
                
        ## We now compute the diameter of the aperture with 50% of the energy
        x=(2*self.rRad+1)*self.pix_scale
        y=profiles[1,]
        t=numpy.greater(y,0.5)
        try:
            i1=numpy.nonzero(t)[0][0]
        except:
            print "ERROR3 in science.py - no non-zero elements... Oh well!"
            i1=0
            
        x1=x[i1]
        t=numpy.less(y,0.5)
        try:
            i2=numpy.nonzero(t)[0][-1]
        except:
            print "ERROR4 in science.py - no non-zero elements... Oh well!"
            i2=0
            
        x2=x[i2]
        if x2==x1:
            x2+=1
        ##we compute the equation of the line between those two points
        a=(y[i2]-y[i1])/(x2-x1)
        b=y[i2]-a*x2
        ##we then compute the diameter
        d50=(0.5-b)/a
        ##print "d50=%g"%(d50)
        dictScience['d50']=d50
        
        ## We now compute the energy into a square aperture of 0.2 arcsec
        ## We compute the ensquared energy profile
        inbox=self.computeEnsquaredEnergy(longExpPSF)
        y=inbox
        x=self.wid
        #apSize=0.2
        for apSize in self.inboxDiamList:
            t=numpy.greater(x,apSize)
            try:
                i1=numpy.nonzero(t)[0][0]
            except:
                print "ERROR5 in science.py - no non-zero elements... Oh well!"
                i1=0
            x1=x[i1]
            t=numpy.less(x,apSize)
            try:
                i2=numpy.nonzero(t)[0][-1]
            except:
                print "ERROR6 in science.py - no non-zero elements... Oh well! inbox %g"%apSize
                i1=0            
            x2=x[i2]
            if x2==x1:
                x2+=1
            ##we compute the equation of the line between those two points
            a=(y[i2]-y[i1])/(x2-x1)
            b=y[i2]-a*x2
            ##we then compute the energy
            inbox=a*apSize+b
            ##print "E(0.2 arcsec)=%g"%(inbox)
            dictScience['inbox%g'%apSize]=inbox

        #now compute speckle fraction.  This is really only meaningful for high strehl cases.
        if self.keepDiffPsf and self.nimg==self.nfft and self.nfft==2*self.npup:
            #find difference
            diff=self.diffPSF-longExpPSF
            #zero out the middle part (up to first null)
            d=diff.shape[0]/2#find the middle
            diff[d-1:d+2,d-1:d+2]=0#and get rid of 13 pixels around it.
            diff[d,d-2:d+3]=0
            diff[d-2:d+3,d]=0
            #and compute the squared sum.
            dictScience['speckle']=(diff**2).sum()

        return dictScience

##     def computeScientificParameters(self):
##         """computeScientificParameters function : computes FWHM, 
##            Modifications made by FA"""
##         ## We do the average of the PSF and normalise it to 1
##         self.longExpPSF[:,]=self.longExpImg/Numeric.sum(Numeric.sum(self.longExpImg))##We normalise the instantaneous PSF to 1
##         #window(2,wait=1)
##         #fma()
##         #pli(self.longExpPSF)
                
##         ## We calculate PSF parameters ###
##         ## We start by the Strehl Ratio
##         nfft=self.nfft
##         strehl=self.longExpPSF[nfft/2,nfft/2]/self.diffn_core_en
##         self.dictScience['strehl']=strehl
##         ##print "Strehl=%g"%(strehl)
##         self.dictScience['strehlPeak']=max(self.longExpPSF.flat)/self.diffn_core_en
##         ## We now compute the FWHM
##         ## First we compute the radial profile
##         ##Pbs with interp function : we have to compute it manually :-((
##         profiles=self.computeRadialProfileAndEncircledEnergy(self.longExpPSF)
##         x=self.rad[:,]
##         y=profiles[0,:,]
##         ##we look for the index of x just before and after y[0]/2
##         t=Numeric.greater(y,y[0]/2)
##         try:
##             i1=Numeric.nonzero(t)[-1]
##         except:
##             print "ERROR in util/science.py - no non-zero elements... Oh well!"
##             i1=0
##         x1=x[i1]
##         t=Numeric.less(y,y[0]/2)
##         try:
##             i2=Numeric.nonzero(t)[0]
##         except:
##             print "ERROR2 in util/science.py - no non-zero elements... Oh well! "
##             i2=0
##         x2=x[i2]
##         ##we compute the equation of the line between those two points
##         if x2==x1:
##             x2+=1
##         a=(y[i2]-y[i1])/(x2-x1)
##         b=y[i2]-a*x2
##         ##we then compute the FWHM
##         fwhm=2*(y[0]/2-b)/a
##         ##print "FWHM=%g"%(fwhm)
##         self.dictScience['FWHM']=fwhm
                
##         ## We now compute the diameter of the aperture with 50% of the energy
##         x=(2*self.rRad+1)*self.pix_scale
##         y=profiles[1,]
##         t=Numeric.greater(y,0.5)
##         try:
##             i1=Numeric.nonzero(t)[0]
##         except:
##             print "ERROR3 in science.py - no non-zero elements... Oh well!"
##             i1=0
            
##         x1=x[i1]
##         t=Numeric.less(y,0.5)
##         try:
##             i2=Numeric.nonzero(t)[-1]
##         except:
##             print "ERROR4 in science.py - no non-zero elements... Oh well!"
##             i2=0
            
##         x2=x[i2]
##         if x2==x1:
##             x2+=1
##         ##we compute the equation of the line between those two points
##         a=(y[i2]-y[i1])/(x2-x1)
##         b=y[i2]-a*x2
##         ##we then compute the diameter
##         d50=(0.5-b)/a
##         ##print "d50=%g"%(d50)
##         self.dictScience['d50']=d50
        
##         ## We now compute the energy into a square aperture of 0.2 arcsec
##         ## We compute the ensquared energy profile
##         self.inbox=self.computeEnsquaredEnergy(self.longExpPSF)
##         y=self.inbox
##         x=self.wid
##         #apSize=0.2
##         for apSize in self.inboxDiamList:
##             t=Numeric.greater(x,apSize)
##             try:
##                 i1=Numeric.nonzero(t)[0]
##             except:
##                 print "ERROR5 in science.py - no non-zero elements... Oh well!"
##                 i1=0
##             x1=x[i1]
##             t=Numeric.less(x,apSize)
##             try:
##                 i2=Numeric.nonzero(t)[-1]
##             except:
##                 print "ERROR6 in science.py - no non-zero elements... Oh well!"
##                 i1=0            
##             x2=x[i2]
##             if x2==x1:
##                 x2+=1
##             ##we compute the equation of the line between those two points
##             a=(y[i2]-y[i1])/(x2-x1)
##             b=y[i2]-a*x2
##             ##we then compute the energy
##             inbox=a*apSize+b
##             ##print "E(0.2 arcsec)=%g"%(inbox)
##             self.dictScience['inbox%g'%apSize]=inbox
##         return self.dictScience


    def conv(self,img):
        """Convolve image with a PSF
        Convolution by an instrumental PSF. Requires to set img
        """
        temp=numpy.fft.rfft2(img)
        convimg=numpy.fft.irfft2(temp*self.sci_psf_fft)
        return convimg
        

##     def conv(self,img):
##         """Convolve image with a PSF
##         Convolution by an instrumental PSF. Requires to set img
##         """
##         temp=FFT.real_fft2d(img)
##         convimg=FFT.inverse_real_fft2d(temp*self.sci_psf_fft)
##         return convimg



    def doScienceCalc(self,inputData,control,curtime=0):
        """compute the science calculation.  Here, inputData is the phase, control is a dictionary of control commands, such as useFPGA, calcRMSFile, zero_science and science_integrate."""
        useFPGA=control["useFPGA"] and self.canUseFPGA
        ##we compute the piston and remove it into the pupil
        if self.atmosPhaseType=="phaseonly":
            if inputData.shape!=(self.npup,self.npup):#are we binning the phase before centroid calculation - might be needed for XAO systems if want to use the fpga (npup max is 1024 for the fpga).
                #This assumes that the inputData size is a power of 2 larger than phs, ie 2, 4, 8 etc times larger.
                #print "Binning pupil for science calculation %d %s "%(self.npup,str(inputData.shape))
                cmod.binimg.binimg(inputData,self.phs)
                self.phs/=self.realPupBinned
                inputData=self.phs#now the binned version.
            pist=numpy.sum(numpy.sum(inputData*self.pup))/self.pupsum ##computation of the piston
            numpy.put(self.phs.ravel(),self.idxPup,numpy.take(inputData.ravel(),self.idxPup)-pist) ##we remove the piston only from stuff in the pupil.
        else:
            raise Exception("science: todo: don't know how to remove piston")
        nfft=self.nfft
        t1=time.time()
        if self.phaseMultiplier!=1:
            self.phs*=self.phaseMultiplier
        if control["calcRMS"]:
            self.phaseRMS=self.calcRMS(self.phs,self.pup)
            self.phaseRMSsum+=self.phaseRMS
            self.phaseRMSsum2+=self.phaseRMS*self.phaseRMS
            self.phaseRMScnt+=1
            if self.calcRMSFile!=None:
                f=open(self.calcRMSFile,"a")
                f.write("%g\t%g\n"%(curtime,self.phaseRMS))
                f.close()
            self.dictScience["rms"]=self.phaseRMS
        else:
            self.dictScience["rms"]=0.

        #if self.phaseTilt!=None:#a phase tilt is needed if binning by an even factor (to put the central spot in a single pixel)
        #    self.phs+=self.phaseTilt
        if control["zero_science"]>0:
            for k in self.dictScience.keys():
                self.dictScience[k]=0.
            self.dictScience["rms"]=self.phaseRMS
            #self.dictScience['FWHM']=0.
            #self.dictScience['d50']=0.
            #self.dictScience['inbox']=0.
            #self.dictScience['strehl']=0.
            self.isamp=0
            self.psfSamp=0
            self.longExpImg[:,]=0.
            self.longExpPSF[:,]=0.
            self.n_integn=0
            self.phaseRMSsum=0.
            self.phaseRMSsum2=0.
            self.phaseRMScnt=0
            self.luckyCnt=0
        if control["science_integrate"]:# Integrate if shutter is open
            if control["zero_science"]==0:# Not zeroing science image
                self.psfSamp+=1
                if self.psfSamp>=self.sciPSFSamp:
                    self.psfSamp=0
                    self.computeShortExposurePSF(self.phs,useFPGA)#calc short exposure PSF
                    self.longExpImg+=self.instImg# Integrate img
                    self.n_integn+=1 ##We increment the number of integrations used to compute long exposure PSF
                    self.isamp+=1##we increment the isamp counter
                    if (self.isamp>=self.scinSamp): #compute scientific parameters
                        self.isamp=0#we reset isamp to 0
                        # We do the average of the PSF and normalise it to 1
                        self.longExpPSF[:,]=self.longExpImg/numpy.sum(self.longExpImg)##We normalise the instantaneous PSF to 1

                        self.computeScientificParameters()
                        # we update the history lists
                        if self.history==None:
                            self.historyKeys=self.dictScience.keys()
                            self.history=numpy.zeros((len(self.historyKeys),self.historyListsSize),numpy.float32)
                            self.historyCnt=0
                        for ik in range(len(self.historyKeys)):
                            k=self.historyKeys[ik]
                            self.history[ik,self.historyCnt]=self.dictScience[k]
                        self.historyCnt=(self.historyCnt+1)%self.history.shape[1]
                    if control["lucky_integrate"]:#doing lucky...
                        if self.luckyCnt==0:
                            if self.luckyImg==None:
                                self.luckyImg=numpy.empty((self.nimg,self.nimg),self.fpDataType)

                            self.luckyImg[:]=self.instImg
                        else:
                            self.luckyImg+=self.instImg
                        self.luckyCnt+=1
                        if self.luckyCnt>=self.luckyNSampFrames:
                            self.luckyCnt=0
                            #Now do the lucky calculations.
                            self.luckyImg/=self.luckyImg.sum()#normalise it
                            self.computeScientificParameters(self.luckyImg,self.luckyDict)
                            if self.luckyHistory==None:
                                self.luckyHistoryKeys=self.luckyDict.keys()
                                self.luckyHistory=numpy.zeros((len(self.luckyHistoryKeys),self.luckyHistorySize),numpy.float32)
                                self.luckyHistoryCnt=0
                            for ik in range(len(self.luckyHistoryKeys)):
                                k=self.luckyHistoryKeys[ik]
                                self.luckyHistory[ik,self.luckyHistoryCnt]=self.luckyDict[k]
                            self.luckyHistoryCnt=(self.luckyHistoryCnt+1)%self.luckyHistory.shape[1]
                            if self.luckyImgFilename!=None and self.luckyImgSize>0:
                                if self.luckyFile==None:
                                    self.luckyFile=tempfile.TemporaryFile()
                                    self.luckyLastImg=numpy.zeros((self.luckyImgSize,self.luckyImgSize),self.luckyImg.dtype)
                                self.luckyLastImg[:]=self.luckyImg[self.nimg/2-self.luckyImgSize/2:self.nimg/2+self.luckyImgSize/2,self.nimg/2-self.luckyImgSize/2:self.nimg/2+self.luckyImgSize/2]
                                self.luckyFile.write(self.luckyLastImg.byteswap().tostring())

                        #if (len(self.strehl_hist)==self.scienceListsSize):
                        #    vout=self.strehl_hist.pop(0)#remove from start of list
                        #    vout=self.fwhm_hist.pop(0)
                        #    vout=self.d50_hist.pop(0)
                        #    vout=self.inbox_hist.pop(0)
                        #self.strehl_hist.append(self.dictScience['strehl'])
                        #self.fwhm_hist.append(self.dictScience['FWHM'])
                        #self.d50_hist.append(self.dictScience['d50'])
                        #if len(self.inboxDiamList)>0:
                        #    self.inbox_hist.append(self.dictScience['inbox%g'%self.inboxDiamList[0]])
##                     if self.fitsFilename!=None:
##                         if os.path.exists(self.fitsFilename):
##                             os.remove(self.fitsFilename);
##                         head=[]
##                         head.append("SCALE   = %g"%self.pix_scale)
##                         head.append("STREHL  = %g"%self.dictScience["strehl"])
##                         head.append("FWHM    = %g"%self.dictScience["FWHM"])
##                         head.append("D50     = %g"%self.dictScience["d50"])
##                         util.FITS.Write(self.longExpImg/self.n_integn,self.fitsFilename,extraHeader=head)
        if(self.timing):print "science",time.time()-t1
        if self.debug!=None:
            print "science: generateNext done (debug=%s)"%str(self.debug)
        return None
##     def doScienceCalc(self,inputData,control,curtime=0):
##         """compute the science calculation.  Here, inputData is the phase, control is a dictionary of control commands, such as useFPGA, calcRMSFile, zero_science and science_integrate."""
##         useFPGA=control["useFPGA"] and self.canUseFPGA
##         ##we compute the piston and remove it into the pupil
##         if self.atmosPhaseType=="phaseonly":
##             if inputData.shape!=(self.npup,self.npup):#are we binning the phase before centroid calculation - might be needed for XAO systems if want to use the fpga (npup max is 1024 for the fpga).
##                 #This assumes that the inputData size is a power of 2 larger than phs, ie 2, 4, 8 etc times larger.
##                 #print "Binning pupil for science calculation %d %s "%(self.npup,str(inputData.shape))
##                 cmod.binimg.binimg(inputData,self.phs)
##                 self.phs/=self.realPupBinned
##                 inputData=self.phs#now the binned version.
##             inputData=Numeric.array(inputData)
##             pist=Numeric.sum(Numeric.sum(inputData*self.pup))/self.pupsum ##computation of the piston

##             Numeric.put(self.phs.flat,self.idxPup,Numeric.take(inputData.flat,self.idxPup)-pist) ##we remove the piston only from stuff in the pupil.
##         else:
##             raise Exception("science: todo: don't know how to remove piston")
##         nfft=self.nfft
##         t1=time.time()
##         if control["calcRMS"]:
##             self.phaseRMS=self.calcRMS(self.phs,self.pup)
##             if self.calcRMSFile!=None:
##                 f=open(self.calcRMSFile,"a")
##                 f.write("%g\t%g\n"%(curtime,self.phaseRMS))
##                 f.close()
##         if control["zero_science"]>0:
##             self.dictScience['FWHM']=0.
##             self.dictScience['d50']=0.
##             self.dictScience['inbox']=0.
##             self.dictScience['strehl']=0.
##             self.isamp=0
##             self.longExpImg[:,]=0.
##             self.n_integn=0
##         if control["science_integrate"]:# Integrate if shutter is open
##             if control["zero_science"]==0:# Not zeroing science image
##                 self.computeShortExposurePSF(self.phs,useFPGA)#calc short exposure PSF
##                 self.longExpImg+=self.instImg# Integrate img
##                 self.n_integn+=1 ##We increment the number of integrations used to compute long exposure PSF
##                 self.isamp+=1##we increment the isamp counter
##                 if (self.isamp>=self.scinSamp): #compute scientific parameters
##                     self.isamp=0#we reset isamp to 0
##                     self.computeScientificParameters()
##                     # we update the history lists
##                     if (len(self.strehl_hist)==self.scienceListsSize):
##                         vout=self.strehl_hist.pop(0)#remove from start of list
##                         vout=self.fwhm_hist.pop(0)
##                         vout=self.d50_hist.pop(0)
##                         vout=self.inbox_hist.pop(0)
##                     self.strehl_hist.append(self.dictScience['strehl'])
##                     self.fwhm_hist.append(self.dictScience['FWHM'])
##                     self.d50_hist.append(self.dictScience['d50'])
##                     if len(self.inboxDiamList)>0:
##                         self.inbox_hist.append(self.dictScience['inbox%g'%self.inboxDiamList[0]])
## ##                     if self.fitsFilename!=None:
## ##                         if os.path.exists(self.fitsFilename):
## ##                             os.remove(self.fitsFilename);
## ##                         head=[]
## ##                         head.append("SCALE   = %g"%self.pix_scale)
## ##                         head.append("STREHL  = %g"%self.dictScience["strehl"])
## ##                         head.append("FWHM    = %g"%self.dictScience["FWHM"])
## ##                         head.append("D50     = %g"%self.dictScience["d50"])
## ##                         util.FITS.Write(self.longExpImg/self.n_integn,self.fitsFilename,extraHeader=head)
##         if(self.timing):print "science",time.time()-t1
##         if self.debug!=None:
##             print "science: generateNext done (debug=%s)"%str(self.debug)
##         return None
class fpgaSciStateInformation:
    def __init__(self):
        self.npup=None
        self.nfft=None
        self.userPrecision=None
        
def difYorick (x, i = 0) :

   """
   Taken from site-packages/gist/yorick.py.
   dif_(x, i) does the same thing as in Yorick: x(...,dif_,...)
   where dif_ is the ith subscript. (works for up to 5 dimensions).
   Namely, the elements along the ith dimension of x are replaced
   by the differences of adjacent pairs, and the dimension decreases
   by one. Remember that Python subscripts are counted from 0.
   """

   if numpy.isscalar (x) :
      raise Exception("util.sci.dif_: must be called with an array.")
   dims = x.shape
   ndims = len (dims)
   if i < 0 or i > ndims - 1 :
      raise Exception("util.sci.dif_: i <" + `i+1` + "> is out of the range of x's dimensions <" + `ndims` +">.")
   if i == 0 :
      newx = x [1:dims [0]] - x [0:dims [0] - 1]
   elif i == 1 :
      newx = x [:, 1:dims [1]] - x[:, 0:dims [1] - 1]
   elif i == 2 :
      newx = x [:, :, 1:dims [2]] - x [:, :, 0:dims [2] - 1]
   elif i == 3 :
      newx = x [:, :, :, 1:dims [3]] - x [:, :, :, 0:dims [3] - 1]
   elif i == 4 :
      newx = x [:, :, :, :, 1:dims [4]] - x [:, :, :, :, 0:dims [4] - 1]
   return newx
