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
import util.atmos
import util.calcPxlScale

class sciOverview:
    def __init__(self,sciDict):
        self.sciDict=sciDict#A dictionary of sciInfo objects
    def getSciByID(self,idstr):
        sci=self.sciDict[idstr]
        if sci.idstr!=idstr:
            raise Exception("Inconsistency in sciOverview")
        return sci
    def values(self):
        return self.sciDict.values()
class sciInfo(util.atmos.source):
    def __init__(self,idstr,theta,phi,pupil,sourcelam,phslam=None,nsamp=10,zeroPsf=10,psfFilename=None,summaryFilename="results.csv",histFilename=None,integrate=1,calcRMS=0,phaseType="phaseonly",nfft=None,nimg=None,realPupil=None,sciPath=None,dmpupil=None,usedmpup=0,psfSamp=1,luckyObj=None,saveString=None,diffPsfFilename=None,histListSize=None,inboxDiamList=[0.2],userFitsHeader=None,psfEnergyToSave=0.,psfMinSize=10,nimgLongExp=None):#Note, if change this, also update the .copy method.
        super(sciInfo,self).__init__(idstr,theta,phi,alt=-1,sourcelam=sourcelam,phslam=phslam)
        self.nsamp=nsamp
        self.zeroPsf=10
        self.psfFilename=psfFilename
        self.summaryFilename=summaryFilename
        self.histFilename=histFilename
        self.integrate=integrate
        self.calcRMS=calcRMS
        self.phaseType=phaseType
        self.nfft=nfft
        self.nimg=nimg
        self.pupil=pupil
        self.realPupil=realPupil
        self.dmpupil=dmpupil
        self.usedmpup=usedmpup
        self.psfSamp=psfSamp
        self.luckyObj=luckyObj#a luckyInfo object.
        self.saveString=saveString
        self.diffPsfFilename=diffPsfFilename
        self.histListSize=histListSize
        self.inboxDiamList=inboxDiamList
        self.userFitsHeader=userFitsHeader
        self.psfEnergyToSave=psfEnergyToSave
        self.psfMinSize=psfMinSize
        self.nimgLongExp=nimgLongExp
        self.sciPath=sciPath
        if self.sciPath is None:
            self.sciPath=idstr
    def copy(self,idstr):
        """Copy to a new ID-string"""
        s=sciInfo(idstr,self.theta,self.phi,self.pupil,self.sourcelam,self.phslam,self.nsamp,self.zeroPsf,self.psfFilename,self.summaryFilename,self.histFilename,self.integrate,self.calcRMS,self.phaseType,self.nfft,self.nimg,self.realPupil,self.sciPath,self.dmpupil,self.usedmpup,self.psfSamp,self.luckyObj,self.saveString,self.diffPsfFilename,self.histListSize,self.inboxDiamList,self.userFitsHeader,self.psfEnergyToSave,self.psfMinSize,self.nimgLongExp)
        return s

class luckyInfo:
    def __init__(self,filename=None,imgFilename=None,imgSize=None,nSampFrames=1,byteswap=0,integrate=1,histSize=None):
        self.filename=filename
        self.imgFilename=imgFilename
        self.imgSize=imgSize
        self.nSampFrames=nSampFrames
        self.byteswap=byteswap
        self.integrate=integrate
        self.histSize=histSize

class science:
    """A class of utilities related to science imaging.  This can be used
    stand-alone or as part of the aosim/science/science.py module.
    realPup - used if phase is binned before image creation.
    Binning phase can be used to reduce computation time, or eg make it fit into the FPGA... for a 42m telescope, 1024x1024 phase values are fine, binning down from 2048 or 4096 gives virtually no difference to strehl ratio.  This is because at 1024, each phase size is about 4cm, which is less than half r0... at 2048, this is even smaller, so doesn't change anything.
    InboxDiamList is a list of diameters (arcsec) for which the ensquared energy should be calculated.
    phaseMultiplier is used for cases where phase wavelength is different from the wavelength that this science object is at.

    """
    def __init__(self,npup,nfft,pup,nimg=None,tstep=0.005,atmosPhaseType="phaseonly",keepDiffPsf=0,pix_scale=1.,fitsFilename=None,diffPsfFilename=None,scinSamp=1,sciPSFSamp=1,scienceListsSize=128,debug=None,timing=0,allocateMem=1,realPup=None,fpDataType="f",calcRMSFile=None,inboxDiamList=[0.2],sciFilename=None,saveFileString=None,nthreads=1,histFilename=None,phaseMultiplier=1,luckyNSampFrames=1,luckyFilename=None,luckyImgFilename=None,luckyImgSize=None,luckyHistorySize=10,luckyByteswap=0,userFitsHeader=None,psfEnergyToSave=0.,psfMinSize=10,nimgLongExp=None):
        self.fftPlan=None
        self.nfft=nfft
        self.npup=npup
        if nimg is None:
            nimg=nfft
        self.nimg=nimg
        if nimgLongExp==None:
            nimgLongExp=self.nimg
        self.nimgLongExp=nimgLongExp
        if (self.nfft%self.nimg)!=0:
            raise Exception("WARNING: util.sci.py nfft not a multiple of nimg")
        self.phaseTilt=None
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
        self.computeDiffPsfProfiles=0#compute radial and ensquared profiles for the diffraction limited psf.
        self.historyListsSize=int(scienceListsSize)
        self.luckyHistorySize=luckyHistorySize
        self.luckyByteswap=luckyByteswap
        self.inboxDiamList=inboxDiamList
        self.saveFileString=saveFileString
        if userFitsHeader is not None and type(userFitsHeader)!=type([]):
            userFitsHeader=[userFitsHeader]
        self.userFitsHeader=userFitsHeader
        self.psfEnergyToSave=psfEnergyToSave
        self.psfMinSize=psfMinSize
        if type(self.saveFileString)!=type(""):
            self.saveFileString=""
        self.debug=debug
        self.nthreads=nthreads
        self.timing=timing
        if type(realPup)!=type(None):#this is used if the phase is binned before image creation.
            self.realPupBinned=numpy.zeros((self.npup,self.npup),numpy.float32)
            cmod.binimg.binimg(realPup.astype("f"),self.realPupBinned)
            self.realPupBinned[:,]=numpy.where(self.realPupBinned==0,1,self.realPupBinned)
        self.fpDataType="f"#always float32 now (cmod.fft requires it).
        self.integratedImgDataType=fpDataType#floating point type
        self.cpDataType='F'#fpDataType.upper()#complex type
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
        #self.longExpPSF=numpy.zeros((nimg,nimg),self.integratedImgDataType)#was float64
        self.longExpImg=numpy.zeros((nimgLongExp,nimgLongExp),self.integratedImgDataType)# Long exposure image array (cannot be shared)
        self.luckyImg=None
        self.luckyCnt=0
        self.luckyHistoryKeys={}
        self.luckyNSampFrames=luckyNSampFrames#Note, every sciPSFSamp (which defaults to 1, not scinsamp) frames, the lucky image will be integrated, until this has been done luckyNSampFrames, at which point results are calculated and the luckyImg zeroed.
        self.luckyDict={}
        self.luckyHistory=None
        self.luckyFilename=luckyFilename
        self.luckyImgFilename=luckyImgFilename
        self.luckyImgSize=luckyImgSize
        self.luckyFile=None
        if allocateMem:
            self.initMem()
            self.initProfiles()
        #self.longExpPSFView=None
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

        
    def initMem(self,fftTmp=None,pupilAmplitude=None,focusAmplitude=None,tempImg=None,instImg=None,phs=None,binImg=None,longExpPSF=None):
        """Initialise memories needed... here, if memories are passed in,
        these are used (ie when sharing resources), otherwise, new memories
        are allocated."""
        nfft=self.nfft
        npup=self.npup
        nimg=self.nimg
        nimgLongExp=self.nimgLongExp
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
        if focusAmplitude is None:#an in place FFT will be done (slower)
            self.focusAmplitude=self.pupilAmplitude
        elif focusAmplitude.shape!=(nfft,nfft):
            self.focusAmplitude=util.arrayFromArray.arrayFromArray(focusAmplitude,(nfft,nfft),self.cpDataType)#Numeric.Complex64
        else:
            self.focusAmplitude=focusAmplitude

        if longExpPSF is None:
            self.longExpPSF=numpy.zeros((nimgLongExp,nimgLongExp),self.integratedImgDataType)
        elif longExpPSF.shape!=(nimgLongExp,nimgLongExp):
            self.longExpPSF=longExpPSF.view(self.integratedImgDataType)[:self.nimgLongExp,:self.nimgLongExp]
        else:
            self.longExpPSF=longExpPSF

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
        
        if type(instImg)==type(None):
            self.instImg=numpy.zeros((self.nimg,self.nimg),self.fpDataType)#was Numeric.Float64
        else:
            #instImg=numpy.array(instImg)#xxx
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

        
    def initProfiles(self):
        """This should be called after the memory has been set up..."""
        self.diffPSF=self.computeDiffPSF()
        self.diffn_core_en=float(self.diffPSF[self.nimg/2,self.nimg/2])
        if self.computeDiffPsfProfiles:
            self.diffPsfRadialProfile=self.computeRadialProfileAndEncircledEnergy(self.diffPSF)[0,:,]
            print "todo - return from computeEnsquaredEnergy - allocate once"
            self.diffPsfEnsquaredProfile=self.computeEnsquaredEnergy(self.diffPSF)
        if self.diffPsfFilename!=None:#save the diffraction limited PSF.
            util.FITS.Write(self.diffPSF,self.diffPsfFilename)
        #self.dlPsfRadialProfile=self.computeRadialProfileAndEncircledEnergy(self.diffPSF)[0,:,]
        if self.computeOTF:
            self.diffOTFSum=numpy.fft.fft2(numpy.fft.fftshift(self.diffPSF),s=(self.diffPSF.shape[0]*2,self.diffPSF.shape[1]*2)).sum()
        #print "Diffraction OTF sum: %s"%str(self.diffOTFSum)
        #perform an acml FFT initialisation.
        #Actually not needed for inplace_fft2d...

    def initRadialProfile(self):
        """Creates and initialises arrays to compute the radial
        profile and the encircled energy profile of the PSF Must be
        first called before calling
        computeRadialProfileAndEncircledEnergy
        """	
        ####  We first create a map of the square of the distances to the center of the PSF image
        ####  The map of distances is computed with the dist function (see aosim/utils/dist.py for help)
        ####  We use the square because we have then integer indexes
        r2=dist(self.nimgLongExp,natype=self.fpDataType,sqrt=0).ravel()
        #r2*=r2 ##no longer sqrts in dist... square : we have integer numbers
        self.nnRad=numpy.argsort(r2) ##we flatten the grid of distances and sort the values
        r2r=numpy.take(r2,self.nnRad) ## we sort the grid of distances
        #self.xRad=numpy.nonzero(difYorick(r2r))[0] ##we look when the grid of distances change of value
        self.xRad=numpy.nonzero(r2r[1:]-r2r[:-1])[0]
        #print self.xRad,type(self.xRad),type(self.xRad[0])
        #self.tabNbPointsRad=difYorick(self.xRad) ## number of points per bin
        self.tabNbPointsRad=self.xRad[1:]-self.xRad[:-1]
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
        #tabEnergyPerPixBin=difYorick(tabEncircledEnergy) ##to have the energy per pixel bin
        tabEnergyPerPixBin=tabEncircledEnergy[1:]-tabEncircledEnergy[:-1]
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


    def initEnsquaredEnergyProfile(self):
        """Creates and initialises arrays to fastly compute the ensquared energy profile
        Must be called before calling computeEnsquaredEnergy
        """	
        ####  We first create a pixel map of concentric square apertures
        tabx=numpy.arange(self.nimgLongExp)-self.nimgLongExp/2;
        r2=numpy.maximum(numpy.absolute(tabx[:,numpy.newaxis]),numpy.absolute(tabx[numpy.newaxis,:,]))

        ##we flatten the grid of distances and sort the values
        self.nnSquare=numpy.argsort(r2.ravel())
        r2r=numpy.take(r2.ravel(),self.nnSquare) ## we sort the grid of distances
        #self.xSquare=numpy.nonzero(difYorick(r2r))[0] ##we look when the grid of distances change of value
        self.xSquare=numpy.nonzero(r2r[1:]-r2r[:-1])[0]
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

    def computeDiffPSF(self):
        """computeDiffPSF function : computes the diffraction limited PSF
        Function written  by FA, heavily modified by AGB"""
        if self.atmosPhaseType=="phaseamp":
            self.computeShortExposurePSF(numpy.array([1,1],self.fpDataType))#places result into instImg
        else:
            self.computeShortExposurePSF(numpy.array([1],self.fpDataType))
        if self.keepDiffPsf:
            self.diffPSF=self.instImg.copy()
        else:
            self.diffPSF=self.instImg#will be overwritten next time computeShortExposurePSF is called...
        return self.diffPSF


    def computeShortExposurePSF(self,phs):
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
                

    def computeScientificParameters(self,longExpPSF=None,dictScience=None):
        """computeScientificParameters function : computes FWHM, 
           Modifications made by FA
           The longExpPSF should be normalised to 1, i.e.sum()==1.
           """
        if longExpPSF is None:
            longExpPSF=self.longExpPSF
        if dictScience is None:
            dictScience=self.dictScience
        #window(2,wait=1)
        #fma()
        #pli(self.longExpPSF)
                
        ## We calculate PSF parameters ###
        ## We start by the Strehl Ratio
        nfft=self.nfft
        nimgLongExp=self.nimgLongExp
        strehl=longExpPSF[nimgLongExp/2,nimgLongExp/2]/self.diffn_core_en
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
        if self.keepDiffPsf and self.nimgLongExp==self.nfft and self.nfft==2*self.npup:
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



    def conv(self,img):
        """Convolve image with a PSF
        Convolution by an instrumental PSF. Requires to set img
        """
        temp=numpy.fft.rfft2(img)
        convimg=numpy.fft.irfft2(temp*self.sci_psf_fft)
        return convimg
        

    def prepareInput(self,inputData):
        """Optionally bins phase and scales for wavelength"""
        ##we compute the piston and remove it into the pupil
        if self.atmosPhaseType=="phaseonly":
            if inputData.shape!=(self.npup,self.npup):#are we binning the phase before centroid calculation - might be needed for XAO systems if want to use the fpga (npup max is 1024 for the fpga).
                #This assumes that the inputData size is a power of 2 larger than phs, ie 2, 4, 8 etc times larger.
                #print "Binning pupil for science calculation %d %s "%(self.npup,str(inputData.shape))
                cmod.binimg.binimg(inputData,self.phs)
                self.phs/=self.realPupBinned
                inputData=self.phs#now the binned version.
            #Don't remove piston any more (160614)
            #pist=numpy.sum(numpy.sum(inputData*self.pup))/self.pupsum ##computation of the piston
            #numpy.put(self.phs.ravel(),self.idxPup,numpy.take(inputData.ravel(),self.idxPup)-pist) ##we remove the piston only from stuff in the pupil.
            else:
                self.phs[:]=inputData
        else:
            raise Exception("science: todo: don't know how to remove piston")
        if self.phaseMultiplier!=1:
            self.phs*=self.phaseMultiplier
        
        
    def doScienceCalc(self,inputData,control,curtime=0):
        """compute the science calculation.  Here, inputData is the phase, control is a dictionary of control commands, such as useFPGA, calcRMSFile, zero_science and science_integrate."""
        ##we compute the piston and remove it into the pupil
        # if self.atmosPhaseType=="phaseonly":
        #     if inputData.shape!=(self.npup,self.npup):#are we binning the phase before centroid calculation - might be needed for XAO systems if want to use the fpga (npup max is 1024 for the fpga).
        #         #This assumes that the inputData size is a power of 2 larger than phs, ie 2, 4, 8 etc times larger.
        #         #print "Binning pupil for science calculation %d %s "%(self.npup,str(inputData.shape))
        #         cmod.binimg.binimg(inputData,self.phs)
        #         self.phs/=self.realPupBinned
        #         inputData=self.phs#now the binned version.
        #     pist=numpy.sum(numpy.sum(inputData*self.pup))/self.pupsum ##computation of the piston
        #     numpy.put(self.phs.ravel(),self.idxPup,numpy.take(inputData.ravel(),self.idxPup)-pist) ##we remove the piston only from stuff in the pupil.
        # else:
        #     raise Exception("science: todo: don't know how to remove piston")
        nfft=self.nfft
        t1=time.time()
        #if self.phaseMultiplier!=1:
        #    self.phs*=self.phaseMultiplier
        if control["calcRMS"]:
            self.prepareInput(inputData)
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
            self.clippedEnergy=0.
            self.n_integn=0
            self.phaseRMSsum=0.
            self.phaseRMSsum2=0.
            self.phaseRMScnt=0
            self.luckyCnt=0
        if control["science_integrate"]:# Integrate if shutter is open
            if control["zero_science"]==0:# Not zeroing science image
                self.psfSamp+=1
                if self.psfSamp>=self.sciPSFSamp:
                    if not control["calcRMS"]:
                        self.prepareInput(inputData)
                    self.psfSamp=0
                    self.computeShortExposurePSF(self.phs)#calc short exposure PSF
                    if self.nimgLongExp==self.nimg:
                        self.longExpImg+=self.instImg# Integrate img
                    else:
                        f=(self.nimg-self.nimgLongExp)//2
                        t=f+self.nimgLongExp
                        self.longExpImg+=self.instImg[f:t,f:t]
                        self.clippedEnergy+=self.instImg[:f].sum()
                        self.clippedEnergy+=self.instImg[t:].sum()
                        self.clippedEnergy+=self.instImg[f:t,:f].sum()
                        self.clippedEnergy+=self.instImg[f:t,t:].sum()
                    #instantaneous strehl calc...
                    self.dictScience['strehlInst']=numpy.max(self.instImg)/self.diffn_core_en
                    self.n_integn+=1 ##We increment the number of integrations used to compute long exposure PSF
                    self.isamp+=1##we increment the isamp counter
                    if (self.isamp>=self.scinSamp): #compute scientific parameters
                        self.isamp=0#we reset isamp to 0
                        # We do the average of the PSF and normalise it to 1
                        self.longExpPSF[:,]=self.longExpImg/(numpy.sum(self.longExpImg)+self.clippedEnergy)##We normalise the instantaneous PSF to 1

                        self.computeScientificParameters()

                        # we update the history lists
                        if self.history is None:
                            self.historyKeys=self.dictScience.keys()
                            self.history=numpy.zeros((len(self.historyKeys),self.historyListsSize),numpy.float32)
                            self.historyCnt=0
                        for ik in range(len(self.historyKeys)):
                            k=self.historyKeys[ik]
                            self.history[ik,self.historyCnt]=self.dictScience[k]
                        self.historyCnt=(self.historyCnt+1)%self.history.shape[1]
                    if control["lucky_integrate"]:#doing lucky...
                        if self.luckyCnt==0:
                            if self.luckyImg is None:
                                self.luckyImg=numpy.empty((self.nimg,self.nimg),self.fpDataType)
                            self.luckyRms=self.phaseRMS
                            self.luckyImg[:]=self.instImg
                        else:
                            self.luckyRms+=self.phaseRMS
                            self.luckyImg+=self.instImg
                        self.luckyCnt+=1
                        if self.luckyCnt>=self.luckyNSampFrames:
                            self.luckyCnt=0
                            #Now do the lucky calculations.
                            self.luckyImg/=self.luckyImg.sum()#normalise it
                            self.computeScientificParameters(self.luckyImg,self.luckyDict)
                            self.luckyDict["rms"]=self.luckyRms/self.luckyNSampFrames#the mean RMS phase that went into this image.
                            if self.luckyHistory is None:
                                self.luckyHistoryKeys=self.luckyDict.keys()
                                self.luckyHistory=numpy.zeros((len(self.luckyHistoryKeys),self.luckyHistorySize),numpy.float32)
                                self.luckyHistoryCnt=0
                            for ik in range(len(self.luckyHistoryKeys)):
                                k=self.luckyHistoryKeys[ik]
                                self.luckyHistory[ik,self.luckyHistoryCnt%self.luckyHistory.shape[1]]=self.luckyDict[k]
                            self.luckyHistoryCnt+=1
                            if self.luckyImgFilename!=None and self.luckyImgSize>0:
                                if self.luckyFile is None:
                                    self.luckyFile=tempfile.TemporaryFile()
                                    self.luckyLastImg=numpy.zeros((self.luckyImgSize,self.luckyImgSize),self.luckyImg.dtype)
                                self.luckyLastImg[:]=self.luckyImg[self.nimg/2-self.luckyImgSize/2:self.nimg/2+self.luckyImgSize/2,self.nimg/2-self.luckyImgSize/2:self.nimg/2+self.luckyImgSize/2]
                                if self.luckyByteswap:
                                    self.luckyLastImg.byteswap(True)
                                self.luckyFile.write(self.luckyLastImg.tostring())
                                if self.luckyByteswap:
                                    self.luckyLastImg.byteswap(True)


        if(self.timing):print "science",time.time()-t1
        if self.debug!=None:
            print "science: generateNext done (debug=%s)"%str(self.debug)
        return None

        
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



def computeScientificParameters(img,nfft=None,nimg=None,npup=None,pupil=None,inboxDiamList=[0.075],lam=1650.,telDiam=4.2,telSec=0.):
    """computeScientificParameters function : computes FWHM, 
       Modifications made by FA
       """
    dictScience={}
    #normalise image...
    img=img/img.sum()


    ## We calculate PSF parameters ###
    ## We start by the Strehl Ratio
    if nfft is None:
        nfft=img.shape[0]
    if nimg is None:
        nimg=img.shape[0]
    if npup is None:
        npup=nfft/2
    if pupil is None:
        pupil=util.tel.Pupil(npup,npup/2,npup/2.*telSec/telDiam)

    #compute pixel scale.
    pix_scale=util.calcPxlScale.pxlScale(lam,telDiam,npup,nfft,nfft/nimg)

    diffPSF=computeShortExposurePSF(0,pupil)
    if diffPSF.shape[0]!=nimg:
        raise Exception("Diffraction limited shape unexpected")
    strehl=img[nimg/2,nimg/2]/diffPSF[nimg/2,nimg/2]
    dictScience['strehl']=strehl
    ##print "Strehl=%g"%(strehl)
    dictScience['strehlPeak']=numpy.max(img)/diffPSF[nimg/2,nimg/2]
    sl=numpy.fft.fftshift(img)
    diffOTFSum=numpy.fft.fft2(numpy.fft.fftshift(diffPSF),s=(diffPSF.shape[0]*2,diffPSF.shape[1]*2)).sum()
    dictScience['strehlOTF']=numpy.fft.fft2(sl,s=(sl.shape[0]*2,sl.shape[1]*2)).sum()/diffOTFSum

    ## We now compute the FWHM
    ## First we compute the radial profile
    ##Pbs with interp function : we have to compute it manually :-((
    rad,rRad,tabNbPointsRad,xRad,nnRad=initRadialProfile(nimg,pix_scale)
    profiles=computeRadialProfileAndEncircledEnergy(img,nnRad,xRad,tabNbPointsRad)
    x=rad
    y=profiles[0,:,]
    ##we look for the index of x just before and after y[0]/2
    t=numpy.greater(y,y[0]/2)
    try:
        i1=numpy.nonzero(t)[0][-1]
    except:
        print "ERROR in util/science.py - no non-zero elements... Oh well! (fwhm invalid)"
        i1=0
    x1=x[i1]
    t=numpy.less(y,y[0]/2)
    try:
        i2=numpy.nonzero(t)[0][0]
    except:
        print "ERROR2 in util/science.py - no non-zero elements... Oh well! (fwhm invalid)"
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
    x=(2*rRad+1)*pix_scale
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
    rSquare,xSquare,nnSquare=initEnsquaredEnergyProfile(nimg)
    wid=rSquare*pix_scale #Square aperture size in arcsec. X Axis of inbox array, used for plotting		
    inbox=computeEnsquaredEnergy(img,nnSquare,xSquare)
    y=inbox

    x=wid
    #apSize=0.2
    for apSize in inboxDiamList:
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


    return dictScience

def computeShortExposurePSF(phs,pup,pad=2):
    """computeShortExposurePSF function : computes short exposure AO corrected PSF
    Modifications made by FA.
    
    If want this to be centred on 4 pixels, rather than 2, add a tilt function to the phs, equal to:
    tiltfn=-(numpy.mgrid[:npup,:npup].sum(0)+1-npup)*numpy.pi/fftsize

    To have a similar function without requiring fliparray2, can add 
    tiltfn*(fftsize+1)
    to the phase (but how does this handle vignetted subaps)
    """

    npup=pup.shape[0] ##to have the dimensions of the input phase array
    pupilAmplitude=numpy.zeros((int(numpy.round(pup.shape[0]*pad)),int(numpy.round(pup.shape[1]*pad))),numpy.complex64)
    pupilAmplitude.real[:npup,:npup]=pup*numpy.cos(phs)
    pupilAmplitude.imag[:npup,:npup]=pup*numpy.sin(phs)
    focAmp=numpy.fft.fft2(pupilAmplitude)
    #img=numpy.absolute(focAmp)
    #img*=img
    #img=focAmp*focAmp.conjugate()
    img=focAmp.real*focAmp.real+focAmp.imag*focAmp.imag
    img=fliparray2(img)
    img/=img.sum()
    return img

def initRadialProfile(nimg,pix_scale):
    """Creates and initialises arrays to compute the radial
    profile and the encircled energy profile of the PSF Must be
    first called before calling
    computeRadialProfileAndEncircledEnergy
    """	
    ####  We first create a map of the square of the distances to the center of the PSF image
    ####  The map of distances is computed with the dist function (see aosim/utils/dist.py for help)
    ####  We use the square because we have then integer indexes
    r2=dist(nimg,natype=numpy.float32,sqrt=0).ravel()
    #r2*=r2 ##no longer sqrts in dist... square : we have integer numbers
    nnRad=numpy.argsort(r2) ##we flatten the grid of distances and sort the values
    r2r=numpy.take(r2,nnRad) ## we sort the grid of distances
    xRad=numpy.nonzero(difYorick(r2r))[0] ##we look when the grid of distances change of value
    #print self.xRad,type(self.xRad),type(self.xRad[0])
    tabNbPointsRad=difYorick(xRad) ## number of points per bin
    rRad=numpy.take(numpy.sqrt(r2r),xRad) ##radius (in pixels) giving the horizontal axis for radial and encircled energy profiles
    rad=rRad*pix_scale
    return rad,rRad,tabNbPointsRad,xRad,nnRad

def initEnsquaredEnergyProfile(nimg):
    """Creates and initialises arrays to fastly compute the ensquared energy profile
    Must be called before calling computeEnsquaredEnergy
    """	
    ####  We first create a pixel map of concentric square apertures
    tabx=numpy.arange(nimg)-nimg/2;
    r2=numpy.maximum(numpy.absolute(tabx[:,numpy.newaxis]),numpy.absolute(tabx[numpy.newaxis,:,]))

    ##we flatten the grid of distances and sort the values
    nnSquare=numpy.argsort(r2.ravel())
    r2r=numpy.take(r2.ravel(),nnSquare) ## we sort the grid of distances
    xSquare=numpy.nonzero(difYorick(r2r))[0] ##we look when the grid of distances change of value
    rSquare=numpy.take(r2r,xSquare)*2+1 ##aperture size (in pixels) giving the horizontal axis for the ensquared energy profile
    return rSquare,xSquare,nnSquare

def computeEnsquaredEnergy(psf,nnSquare,xSquare):
    """Computes ensquared energy profile
    """
    #### We flatten the diffraction limited PSF and sort the values
    psfSort=numpy.take(psf.ravel(),nnSquare)
    #### We compute the encircled energy of the diffraction limited PSF 
    tabEnsquaredEnergy=numpy.take(numpy.cumsum(psfSort),xSquare)
    #### We convert the profile into psf type and return the array
    result=tabEnsquaredEnergy.astype(psf.dtype)
    return result

def computeRadialProfileAndEncircledEnergy(psf,nnRad,xRad,tabNbPointsRad):
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
    psfSort=numpy.take(psf.ravel(),nnRad)
    #### We compute the encircled energy of the PSF
    tabEncircledEnergy=numpy.take(numpy.cumsum(psfSort),xRad)
    #### We compute the radial profile
    tabEnergyPerPixBin=difYorick(tabEncircledEnergy) ##to have the energy per pixel bin
    #tabEnergyPerPixBin.savespace(1)#prevent conversion to double.
    profil=tabEnergyPerPixBin/tabNbPointsRad

    #### Allocation of the return result
    result=numpy.zeros((2,len(tabEncircledEnergy)),dtype=psf.dtype)#savespace=1)

    #### We first store the radial profile in line 1 of the returned array
    result[0,0]=psfSort[0]
    result[0,1:,]=profil[:,] ##radial profile of PSF

    #### We  store the encircled energy profile in line 2 of the returned array
    result[1,:,]=tabEncircledEnergy[:,]
    return result
