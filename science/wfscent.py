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
"""A module for simulating a WFS and centroider, this can make use of the FPGAs, offering speedups of up to 400 times."""
#$Id: wfscent.py,v 1.65 2011/12/06 05:14:39 ali Exp $
### Shack Hartmann WFS simulation object  #####
# Integrate WFS CCD images

import math
import numpy,numpy.random
#import cmod.imgnoise
#from imgnoise import *
import util.centroid
import util.guideStar
#import cmod.binimg
#import cmod.mkimg
#import cmod.mkimgfloat
#import cmod.utils
#import thread
import time,string
import base.aobase
## haveFPGA=1
## try:
##     import fpga
## except:
##     print "FPGA module not installed."
##     haveFPGA=0

class wfscent(base.aobase.aobase):
    """WFS and centroider simulation object.
    This class can make use of resource sharing (memory and FPGA use), by passing a list
    of parents and idstr.
    Shared memory includes:
    self.subimg (hi light images after FFTs)
    self.bimg (SHS images - binned with noise added)
    self.shimg (for display)  - also now centDisplay
    self.pupsub - if pupil.fn is the same (partially shared).
    outputData
    NOT self.reorderedPhs - unless only 1 integration.
    Not sure whether the cell stuff can be used with resource sharing - probably not, since reorderedPhs can't be shared (unless only 1 integ).
    Note, not a good idea to use resource sharing, unless you make a copy of outputData as soon as it is produced...

    Can take 1 or 2 parents, if 2, one is selected depending on the value
    of control["parentType"], which can be "open" or "closed".  This is
    useful for switching between open and closed loop operation, eg initially
    closed for poking, then open for running - the parent in this case would
    be a DM object and an atmos object.

    For pixel scales:
    wfs_nsubx*wfs_n does not need to be equal to npup - it will automatically be interpolated from the npup phase pixels.
    If using a psf, it is probably necessary to zeropad (unless the psf is small).  However, some binning can be done before application of the psf, so that the psf does not need to be twice the fft size.  The basic principle is:
    wfs_n.  Pad up to
    fftsize.  Perform fft.
    Optionally pre-bin.
    pad up to psf size (typically x2).
    convolve.
    clip to clipsize.
    bin to nimg.
    
    Note, psfsize must not be greater than fftsize if prebinning.
    
    """

    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Initialise the WFS simulation object
        Parent - parent object to obtain SH images from
        Config - configuration object
        args - Dictionary of additional arguemnts (passed to aobase).
        Also, looks for "fpprecision", which can be 32 or 64, to specify
        floating point precision.
        
        Expects from config file (list incomplete):
         - npup - n pupils
         - wfs_nsubx - number of SH spots in each direction.
         - timing - flag - whether to print timing messages
         - npup - number of pupils
         - telDiam - Diameter of telescope
         - wfs_nsubx - number of SH spots in each direction
         - wfs_n - number of phase values per subap (in 1D).
         - wfs_nfft - FFT size for SH image creation
         - wfs_nimg - SH image size
         - wfs_ncen - number of image pixels to use (central)
         - wfs_floor - Floor to subtract from CCD images
         - pupil - pupil object
         - wfs_minarea - min unvignetted subap area to use
        
        @param parent: The parent object
        @type parent: Object
        @param config: Configuration object
        @type config: Instance
        @param args: Optional arguments
        @type args: Dict
        """
        if type(parent)!=type({}):
            parent={"closed":parent}
            
        #Note - wfscent can now take parent as a dictionary, with keys
        #either open or closed.  Most simulations will just continue
        #to pass a parent object rather than a dictionary, but open
        #loop simulations, which also are used closed loop to generate
        #a poke matrix will take as a dictionary.  The wfscent can
        #then change its parent depending on control...
            
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        obj=self.config.getVal("wfsOverview",raiseerror=0)
        self.obj=obj
        if obj==None:
            print("DEPRECATION: warning: wfsOverview not specified")
        self.imageOnly=self.config.getVal("imageOnly",default=0)# 0 to return slopes, 
                                                                # 1 to return image as nsubx,nsubx,nimg,nimg, 
                                                                # 2 to return image as a 2d image.
                                                                
        if forGUISetup==1:#this won't be set if resource sharing - whats the point!  So, this is ok.
            fullOutput=self.config.getVal("fullWFSOutput",default=1)
            if obj==None:
                nsubx=self.config.getVal("wfs_nsubx")
                nimg=self.config.getVal("wfs_nimg")
            else:
                wfsobj=obj.getWfsByID(idstr)
                nsubx=wfsobj.nsubx
                nimg=wfsobj.nimg
            if fullOutput==1:
                if self.imageOnly==0:
                    self.outputData=[(nsubx,nsubx,2),numpy.float32]
                elif self.imageOnly==1:
                    self.outputData=[(nsubx,nsubx,nimg,nimg),numpy.float32]
                else:
                    self.outputData=[(nsubx*nimg,nsubx*nimg),numpy.float32]
            else:
                if obj==None:
                    n=self.config.getVal("wfs_n")
                    wfs_minarea=self.config.getVal("wfs_minarea")
                    pupfn=self.config.getVal("pupil")
                else:
                    n=wfsobj.phasesize
                    wfs_minarea=wfsobj.minarea
                    pupfn=wfsobj.pupil
                if type(pupfn)!=numpy.ndarray:
                    pupfn=pupfn.fn
                pupfn=pupfn.astype(numpy.float32)
                frac=wfs_minarea*n*n
                nsubaps=0
                for i in xrange(nsubx):        
                    for j in xrange(nsubx):
                        if pupfn[i*n:(i+1)*n,j*n:(j+1)*n].sum()>frac:
                            nsubaps+=1
                self.outputData=[(nsubaps,2),numpy.float32]
        else:
            self.debug=debug
            useCmod=self.config.getVal("useCmod",default=1)
            self.printLinearisationForcing=self.config.getVal("printLinearisationForcing",default=0)
            self.fullOutput=self.config.getVal("fullWFSOutput",default=1)

            if useCmod==0:#only need a random seed if using the python version
                seed=config.getVal("wfs_seed",default=self.config.getVal("seed",raiseerror=0,default=0),raiseerror=0)
                if seed==None:
                    seed=0
                if seed==0:
                    numpy.random.seed(None)
                else:
                    numpy.random.seed(seed)
            calSource=self.config.getVal("calsource",default=0)
            self.fpDataType=numpy.float32#self.config.getVal("fpDataType",default=numpy.float32)
            self.doneFinalInit=0
            #seed here - no point seeding each resource sharing object...
            #random.seed(config.getVal("wfs_seed",default=0))
            #cmod.imgnoise.seed(config.getVal("wfs_rand_seed",default=self.config.getVal("seed",default=0)))
            #cmod.imgnoise.seed(seed)
            # if args.has_key("fpprecision"):
            #     if args["fpprecision"]==32:
            #         self.fpDataType=numpy.float32
            #     else:
            #         self.fpDataType=numpy.float64
            self.sentPlotsCnt=0
            if not useCmod:
               print("INFORMATION:wfscent: Using **Python** implementation")
            #self.cmodcentseed=seed#self.config.getVal("cmodcentseed",default=self.config.getVal("seed",default=0))
            self.nthreads=self.config.getVal("nthreads",default="all")#usually an integer... or "all"
            if self.nthreads=="all":#use all available CPUs...
                self.nthreads=self.config.getVal("ncpu")#getCpus()
                print("INFORMATION:wfscent: Using {0:d} threads".format(self.nthreads))
##(old)                print("INFORMATION:wfscent: Using %d threads"%self.nthreads)
            self.imgmask=1#value used when creating sh img in drawCents()
            self.control={"cal_source":calSource,"useCmod":useCmod,"zeroOutput":0,"parentType":"closed"}#ZeroOutput can be used when calibrating a pseudo-open loop system.
            self.lastCalSource=0#resourcesharing - todo...?


            self.timing=self.config.getVal("timing",default=0)
            for i in xrange(len(self.idstr)):
                idstr=self.idstr[i]
                parent=self.parentList[i]
                self.initialise(parent,idstr)

    def initialise(self,parent,idstr):
        """Perform initialisation for this resource-sharing object."""
        if type(parent)!=type({}):
            parent={"closed":parent}

        this=base.aobase.resourceSharer(parent,self.config,idstr,self.moduleName)
        this.wfsobj=None
        if self.obj!=None:
            this.wfsobj=self.obj.getWfsByID(idstr)
        wfsobj=this.wfsobj
        print("INFORMATION:wfscent:Initialising centroid object {0:s} {1:s}".format(
              str(parent),str(idstr)) )
##(old)        print("INFORMATION:wfscent:Initialising centroid object %s %s"%(str(parent),str(idstr))
        self.thisObjList.append(this)
        #fpDataType=this.config.getVal("fpDataType",default=Numeric.Float32) # Numeric.Float64
        #if self.fpDataType!=None:#specified by args{}
        #    fpDataType=self.fpDataType
        tstep=this.config.getVal("tstep")                                  # Simulation tick length
        telDiam=this.config.getVal("telDiam")                              # Telescope aperture diameter
        this.laserGuideStar=None

        if wfsobj==None:
            atmosPhaseType=this.config.getVal("atmosPhaseType",default="phaseonly")
            # No. of subaps across tel. pupil
            wfs_nsubx=this.config.getVal("wfs_nsubx")                          
            # Phase array size for subaps (phasesize)
            wfs_n=this.config.getVal("wfs_n")                                  
            # FFT size for subap image calc (fftsize)
            wfs_nfft=this.config.getVal("wfs_nfft",default=wfs_n*2)
            clipsize=this.config.getVal("wfs_clipsize",default=wfs_nfft)
            # Image size for subap (pixels) (ie size after binning).
            wfs_nimg=this.config.getVal("wfs_nimg",default=clipsize/2)
            wfs_int=this.config.getVal("wfs_int",tstep)                              # WFS integration time
            integstepFn=None
            wfs_rowint=this.config.getVal("wfs_rowint",default=None,raiseerror=0)           # row integration time
            wfs_read_mean=this.config.getVal("wfs_read_mean",0.)                  # WFS Readnoise e-

            wfs_read_sigma=this.config.getVal("wfs_read_sigma",0.)
            threshType=this.config.getVal("threshType",default=0)
            wfs_lat=this.config.getVal("wfs_lat",0)                              # WFS readout latency
            skybrightness=this.config.getVal("wfs_skybrightness",0.)
            pupil=this.config.getVal("pupil")
            wfs_minarea=this.config.getVal("wfs_minarea",0.5)#0.5... # Min unvignetted subap area to use - why here ????
            wfs_ncen=this.config.getVal("wfs_ncen",default=wfs_nimg)            # Centroiding box size (pixels)
            wfs_floor=this.config.getVal("wfs_floor",wfs_read_mean+skybrightness)          # Centroiding floor value
            sig=this.config.getVal("wfs_sig",default=None,raiseerror=0)
            if sig==None:
                wfs_bandwidth=this.config.getVal("wfs_bandwidth")                  # WFS Optical bandwidth
                wfs_thruput=this.config.getVal("wfs_thruput")                      # Thruput to the WFS
                wfs_phot_rate_factor=this.config.getVal("wfs_phot_rate_factor")    # Photons/sec/area/angstrom
                wfs_mag=this.config.getVal("wfs_mag")                              # Guide star magnitude
                sig=util.centroid.wfs_sig(wfs_bandwidth,wfs_thruput,wfs_phot_rate_factor,telDiam,wfs_nsubx,wfs_int,wfs_mag)
            addPoisson=this.config.getVal("addPoisson",default=1)
            laserGuideStar=this.config.getVal("laserGuideStar",default=None,raiseerror=0)#this can be an instance of util.guideStar.wfs or util.elong.make(...)[0]
            opticalBinning=this.config.getVal("opticalBinning",default=0)#whether to simulate optical binning (beam splitter and 2 cylindrical lenslet arrays and 2 detectors)
            magicCentroiding=this.config.getVal("magicCentroiding",default=0)
            linearSteps=this.config.getVal("linearSteps",default=None,raiseerror=0)
            calNCoeff=this.config.getVal("calNCoeff",default=0)
            stepRangeFrac=this.config.getVal("stepRangeFrac",default=1.)

            centWeight=this.config.getVal("centWeight",raiseerror=0)
            correlationCentroiding=this.config.getVal("correlationCentroiding",default=0)
            if correlationCentroiding:
                corrThresh=this.config.getVal("corrThresh")
                corrPattern=this.config.getVal("corrPattern",raiseerror=0)
                if type(corrPattern)==type(""):
                    corrPattern=util.FITS.Read(corrPattern)[1]
            else:
                corrThresh=0.
                corrPattern=None
                # imageOnly=this.config.getVal("imageOnly",default=0)#0 to return slopes, 1 to return image as nsubx,nsubx,nimg,nimg, and 2 to return image as a 2d image.
            useBrightest=this.config.getVal("useBrightest",default=0)
            preBinningFactor=this.config.getVal("preBinningFactor",default=1)
            parabolicFit=0#use wfsObj if you want these!
            gaussianFitVals=None
            seed=this.config.getVal("wfs_seed",default=this.config.getVal("seed",default=0,raiseerror=0),raiseerror=0)
            subapLocation=None
            centroidPower=None
        else:
            atmosPhaseType=wfsobj.atmosPhaseType
            wfs_nsubx=wfsobj.nsubx
            wfs_n=wfsobj.phasesize
            wfs_nfft=wfsobj.nfft
            clipsize=wfsobj.clipsize
            wfs_nimg=wfsobj.nimg
            wfs_int=wfsobj.integSteps*tstep
            integstepFn=wfsobj.integstepFn#can be None or a function that returns an integer value not greater than integSteps.
            wfs_rowint=wfsobj.rowint
            if wfs_rowint!=None:
                wfs_rowint*=tstep
            
            wfs_read_mean=wfsobj.bglevel
            wfs_read_sigma=wfsobj.readoutNoise
            threshType=wfsobj.threshType
            wfs_lat=wfsobj.latency*tstep
            skybrightness=wfsobj.skyBrightness
            pupil=wfsobj.pupil
            wfs_minarea=wfsobj.minarea
            wfs_ncen=wfsobj.ncen
            wfs_floor=wfsobj.floor
            sig=wfsobj.sig
            addPoisson=wfsobj.addPoisson
            laserGuideStar=wfsobj.lgsPsf#are we using a LGS?
            opticalBinning=wfsobj.opticalBinning
            magicCentroiding=wfsobj.magicCentroiding
            linearSteps=wfsobj.linearSteps
            calNCoeff=wfsobj.calNCoeff
            stepRangeFrac=wfsobj.stepRangeFrac
            centWeight=wfsobj.centWeight
            correlationCentroiding=wfsobj.correlationCentroiding
            corrThresh=wfsobj.corrThresh
            corrPattern=wfsobj.corrPattern
            useBrightest=wfsobj.useBrightest
            preBinningFactor=wfsobj.preBinningFactor
            parabolicFit=wfsobj.parabolicFit
            gaussianFitVals=wfsobj.gaussianFitVals#None, or a tuple of 2 values: (gaussianMinVal, gaussianReplaceVal).
            seed=wfsobj.seed
            subapLocation=wfsobj.subapLocation
            centroidPower=wfsobj.centroidPower
        this.atmosPhaseType=atmosPhaseType
        if atmosPhaseType not in ["phaseonly","phaseamp","realimag"]:
            raise Exception("wfscent: atmosPhaseType not known %s"%atmosPhaseType)
        #todo: this bit should really be in the centroid module. Not needed here...
        bin=0
        if(clipsize!=wfs_nimg):
            bin=1
            if (float(clipsize)/wfs_nimg)%1!=0:
                print("WARNING:wfscent: Non-integer binning of image - "+
                     "software (binimgmodule.so) cannot handle")
        nIntegrations=int(math.ceil(wfs_int/tstep))        

        if (wfs_int/tstep)%1!=0:
            #print "Warning: wfs - Integration times is not a whole number of timesteps - you might misinterpret the results... %g %g"%(wfs_int,tstep)
            print(("WARNING:wfscent: Integration times is not a whole number "+
                  "of timesteps - you might misinterpret the results... {0:g} "+
                  "{1:g}").format(wfs_int,tstep))
        if wfs_rowint!=None and wfs_rowint%tstep!=0:
            print(("ERROR: wfscent: Row integration times is not a whole "+
                  "number of timesteps - results will be unreliable"+
                  ": wfs_rowint%tstep={0:g}").format(wfs_rowint%tstep))



        
        #Now look to see if we're using LGS...
        if laserGuideStar is None:
            lgsType=this.config.getVal("lgsType",default=None,raiseerror=0)#can be None, "sodium" or "rayleigh".
            if lgsType in ["sodium","rayleigh"]:
                # some extra configuration is needed for setting up the spot elongation...
                lgsPower=this.config.getVal("lgsPower",default=5.)
                telFocalLen=this.config.getVal("telFocalLen",default=telDiam*15.)
                atmosGeom=this.config.getVal("atmosGeom",raiseerror=0)
                if atmosGeom==None:
                    r0=this.config.getVal("r0",default=0.2)
                else:
                    r0=atmosGeom.r0
                subapFocalLen=this.config.getVal("subapFocalLen",default=209.)
                subapFov=this.config.getVal("subapFov",default=14.)
                fluxObj=util.guideStar.Flux(wavelength={"sodium":589.,"rayleigh":514.}[lgsType],power=lgsPower,frame_rate=1./(wfs_int+wfs_lat))
                laserGuideStar=util.guideStar.wfs(pup=pupil.fn,nsubx=wfs_nsubx,subap_pix=clipsize,subap_fov=subapFov,tel_diam=telDiam,r0=r0,tel_fl=telFocalLen,subap_fl=subapFocalLen,fluxObj=fluxObj,lgsDefault=lgsType)
        spotpsf=None
        if laserGuideStar is not None:
            if type(laserGuideStar)==type(""):
                laserGuideStar=util.FITS.Read(laserGuideStar)[1]
            if type(laserGuideStar)==numpy.ndarray:
                if laserGuideStar.shape[:2]!=(wfs_nsubx,wfs_nsubx) or laserGuideStar.shape[2]<clipsize or laserGuideStar.shape[3]<clipsize:
                    raise Exception("Laserguide star image wrong shape %s %d %d"%(str(laserGuideStar.shape),wfs_nsubx,clipsize))
                sig=numpy.sum(numpy.sum(laserGuideStar,3),2).astype(numpy.float32)
                spotpsf=laserGuideStar
                this.laserGuideStar=laserGuideStar
            else:
                if laserGuideStar.subap_pix<clipsize or laserGuideStar.nsubx!=wfs_nsubx:
                    raise Exception("laserGuideStar sizes incorrect (clipsize>%g or nsubx!=%g)"%(laserGuideStar.subap_pix,laserGuideStar.nsubx))
                print("INFORMATION:wfscent:Generating LGS spot patterns... (if these are all the same and take a while, you might be better specifying this in the param file...)")
                laserGuideStar.wfs_image(off=0.5)
                laserGuideStar.sumSubaps()
                print("INFORMATION:wfscent:Finished generating LGS spot patterns.")
                spotpsf=laserGuideStar.spotPsf
                sig=laserGuideStar.sig
                this.laserGuideStar=laserGuideStar
            
        subtractTipTilt=this.config.getVal("subtractTipTilt",default=int(this.laserGuideStar is not None),warn=1)
        if type(sig)==type(0.):
            #print("INFORMATION:wfscent:sig is {0:g} phot/subap".format(sig))
##(old)            print "INFORMATION:wfscent:sig is %g phot/subap"%sig
            pass
        else:
            print("INFORMATION:wfscent:sig is array with max %g phot/subap"%max(sig.flat))
##(old)            print "INFORMATION:wfscent:sig is array with max %g phot/subap"%max(sig.flat)
        if type(spotpsf)==type(None):
            if wfsobj==None:
                spotpsf=this.config.getVal("spotpsf",default=None,raiseerror=0)#a spot PSF, eg an airy disc (eg createAiryDisc(self.fftsize,self.fftsize/2,0.5,0.5)), or LGS elongated spots.
            else:
                spotpsf=wfsobj.spotpsf
        if wfsobj.cameraImage:
            if parent.values()[0].outputData.dtype.char!="f" or parent.values()[0].outputData.flags.contiguous==False:
                raise Exception("wfscent with cameraImage!=0 requires a float32 contiguous input from parent")
            inputImage=parent.values()[0].outputData
        else:
            inputImage=None
            
        phslam=None
        wfslam=None
        atmosGeom=this.config.getVal("atmosGeom",raiseerror=0)
        if atmosGeom!=None:
            if wfsobj==None:
                wfslam=atmosGeom.sourceLambda(idstr)
            else:
                wfslam=wfsobj.sourcelam
            phslam=atmosGeom.phaseLambda(idstr)
        if wfslam==None:
            wfslam=this.config.getVal("sourceLam")
            print("WARNING:wfscent: wfslam not in atmosGeom/wfsOverview, using {0:g}".format(wfslam))
##(old)            print("WARNING:wfscent: wfslam not in atmosGeom, using %g"%wfslam)
        if phslam==None:
            phslam=wfslam
        phaseMultiplier=phslam/wfslam
        if phaseMultiplier!=1:
            print("WARNING:wfscent: {0:s} has phaseMultiplier of {1:g}".format(
                  str(idstr), phaseMultiplier) )
##(old)            print "Warning - wfscent %s has phaseMultiplier of %g"%(str(idstr),phaseMultiplier)
#        this.wfscentObj=util.centroid.centroid(wfs_nsubx,pup=pupil,oversamplefactor=None,readnoise=wfs_read_sigma,readbg=wfs_read_mean,addPoisson=1,noiseFloor=wfs_floor,binfactor=None,sig=sig,skybrightness=skybrightness,warnOverflow=None,atmosPhaseType=atmosPhaseType,fpDataType=self.fpDataType,useFPGA=useFPGA,waitFPGA=waitFPGA,waitFPGATime=waitFPGATime,phasesize=wfs_n,fftsize=wfs_nfft,clipsize=clipsize,nimg=wfs_nimg,ncen=wfs_ncen,tstep=tstep,integtime=wfs_int,latency=wfs_lat,wfs_minarea=wfs_minarea,spotpsf=spotpsf,opticalBinning=opticalBinning,useCell=self.control["useCell"],waitCell=1,usecmod=self.control["useCmod"],subtractTipTilt=subtractTipTilt,magicCentroiding=magicCentroiding,linearSteps=linearSteps,stepRangeFrac=stepRangeFrac,phaseMultiplier=phaseMultiplier,centWeight=centWeight,correlationCentroiding=correlationCentroiding,corrThresh=corrThresh,corrPattern=corrPattern,threshType=threshType,imageOnly=self.imageOnly,calNCoeff=calNCoeff,useBrightest=useBrightest)
        #print(wfs_nsubx,addPoisson,atmosPhaseType,calNCoeff,centWeight,clipsize,correlationCentroiding,corrPattern,corrThresh,wfs_nfft,self.fpDataType,self.imageOnly,wfs_int,wfs_lat,linearSteps,magicCentroiding,wfs_ncen,wfs_nimg,wfs_floor,opticalBinning,phaseMultiplier,wfs_n,self.printLinearisationForcing,pupil,wfs_read_mean,wfs_read_sigma,wfs_rowint,sig,skybrightness,spotpsf,stepRangeFrac,subtractTipTilt,threshType,tstep,useBrightest,self.control["useCmod"],wfs_minarea,preBinningFactor)

        this.wfscentObj=util.centroid.centroid(
            wfs_nsubx,
            addPoisson=addPoisson,
            atmosPhaseType=atmosPhaseType,
            binfactor=None,
            calNCoeff=calNCoeff,
            centWeight=centWeight,
            clipsize=clipsize,
            correlationCentroiding=correlationCentroiding,
            corrPattern=corrPattern,
            corrThresh=corrThresh,
            fftsize=wfs_nfft,
            fpDataType=self.fpDataType,
            imageOnly=self.imageOnly,
            integtime=wfs_int,
            latency=wfs_lat,
            linearSteps=linearSteps,
            magicCentroiding=magicCentroiding,
            ncen=wfs_ncen,
            nimg=wfs_nimg,
            noiseFloor=wfs_floor,
            opticalBinning=opticalBinning,
            oversamplefactor=None,
            phaseMultiplier=phaseMultiplier,
            phasesize=wfs_n,
            printLinearisationForcing=self.printLinearisationForcing,
            pup=pupil,
            readbg=wfs_read_mean,
            readnoise=wfs_read_sigma,
            rowintegtime=wfs_rowint,
            sig=sig,
            skybrightness=skybrightness,
            spotpsf=spotpsf,
            stepRangeFrac=stepRangeFrac,
            subtractTipTilt=subtractTipTilt,
            threshType=threshType,
            tstep=tstep,
            useBrightest=useBrightest,
            usecmod=self.control["useCmod"],
            warnOverflow=None,
            wfs_minarea=wfs_minarea,
            preBinningFactor=preBinningFactor,
            parabolicFit=parabolicFit,
            gaussianFitVals=gaussianFitVals,
            seed=seed,
            integstepFn=integstepFn,
            inputImage=inputImage,
            subapLocation=subapLocation,
            centroidPower=centroidPower,
        )


    def finalInitialisation(self):
        """This gets called just before the main loop is entered - to finalise setting up of arrays.  We now know how large they should be, since all resource shares have been added.
        """
        #first find the maximum sizes for the arrays, and whether any objects
        #can use the FPGAs.
        #print "wfscent finalInitialisation"
        if self.doneFinalInit:
            return
        self.doneFinalInit=1
        fftsizemax=0
        phasesizemax=0
        nsubxmax=0
        nsubxnfftmax=0
        maxarrsize=0
        nIntMax=0
        maxnpup=0
        rpmax=0
        imgsizemax=0
        corrimgsizemax=0
        atmosfactor=1#phaseonly...
        canSharePhs=1
        canSharePupsub=1
        pup=None
        readouttime=None
        for this in self.thisObjList:#find the maximum array size...
            wfs=this.wfscentObj
            so=["wfscent","globals"]
            if this.idstr!=None and len(this.idstr)>0:
                this.config.setSearchOrder(["wfscent_%s"%this.idstr]+so)
            else:
                this.config.setSearchOrder(so)
            if wfs.fftsize>fftsizemax:
                fftsizemax=wfs.fftsize
            if wfs.nsubx>nsubxmax:
                nsubxmax=wfs.nsubx
            if wfs.phasesize>phasesizemax:
                phasesizemax=wfs.phasesize
            if wfs.nsubx*wfs.nimg>imgsizemax:
                imgsizemax=wfs.nsubx*wfs.nimg
            if wfs.corrPattern is not None and wfs.corrPattern.shape[-2]*wfs.nsubx>corrimgsizemax:
                corrimgsizemax=wfs.corrPattern.shape[-2]*wfs.nsubx
            if wfs.nsubx*wfs.fftsize>nsubxnfftmax:
                nsubxnfftmax=wfs.nsubx*wfs.fftsize
            if wfs.nsubx*wfs.phasesize>maxnpup:
                maxnpup=wfs.nsubx*wfs.phasesize
            if wfs.nsubx*wfs.nsubx*wfs.nIntegrations*wfs.phasesize*wfs.phasesize>rpmax:
                atmosPhaseType=this.atmosPhaseType#config.getVal("atmosPhaseType",default="phaseonly")
                #print("INFORMATION:wfscent: this.config.searchOrder='{0:s}'".format(str(this.config.searchOrder)) )
##(old)                print this.config.searchOrder
                if atmosPhaseType=="phaseonly":
                    rpmax=wfs.nsubx*wfs.nsubx*wfs.nIntegrations*wfs.phasesize*wfs.phasesize
                else:
                    rpmax=2*wfs.nsubx*wfs.nsubx*wfs.nIntegrations*wfs.phasesize*wfs.phasesize
            if wfs.nIntegrations>nIntMax:
                nIntMax=wfs.nIntegrations
            if readouttime==None:
                readouttime=wfs.integtime+wfs.latency
            if type(pup)==type(None):
                pup=wfs.pupfn
            if wfs.pupfn is not pup:#is the pupil the same?
                if wfs.pupfn.shape==pup.shape:
                    if numpy.sum((wfs.pupfn==pup).flat)!=reduce(lambda x,y:x*y,pup.shape):
                        canSharePupsub=0
                else:
                    canSharePupsub=0
            if wfs.integtime>wfs.tstep:
                canSharePhs=0
            if wfs.integtime+wfs.latency!=readouttime:
                canSharePhs=0
            arrsize=wfs.nsubx*wfs.nsubx*wfs.nIntegrations*wfs.phasesize*wfs.phasesize*4*atmosfactor+wfs.nsubx*wfs.nsubx*2*4
            if arrsize>maxarrsize:
                maxarrsize=arrsize
        if corrimgsizemax<imgsizemax:
            corrimgsizemax=imgsizemax
        #now allocate memories...
        reorderedPhsMem=None
        if canSharePhs:
            reorderedPhsMem=numpy.zeros((rpmax,),numpy.float32)
        if self.imageOnly:
            outputDataMem=numpy.zeros((imgsizemax*imgsizemax,),numpy.float32)
        else:
            outputDataMem=numpy.zeros((nsubxmax,nsubxmax,2),numpy.float32)
        bimgMem=numpy.zeros((imgsizemax*imgsizemax,),numpy.float64)#needs 64 bits for cmod.binimg...
        self.shimg=numpy.zeros((corrimgsizemax,corrimgsizemax),numpy.float32)
        subimgMem=numpy.zeros((nsubxnfftmax*nsubxnfftmax),numpy.float64) 
        if canSharePupsub:
            pupsubMem=numpy.zeros((maxnpup*maxnpup),self.fpDataType)
        else:
            pupsubMem=None

        #now init each wfscentObj with this mem.
        for this in self.thisObjList:
            wfs=this.wfscentObj
            #set up the memories
            #wfs.initMem(fpgarequired,fpgaarr=fpgaarr,shareReorderedPhs=canSharePhs,reorderedPhsMem=reorderedPhsMem,subimgMem=subimgMem,bimgMem=bimgMem,pupsubMem=pupsubMem,outputDataMem=outputDataMem)
            wfs.initMem(shareReorderedPhs=canSharePhs,reorderedPhsMem=reorderedPhsMem,subimgMem=subimgMem,bimgMem=bimgMem,pupsubMem=pupsubMem,outputDataMem=outputDataMem)
            wfs.finishInit()
            wfs.initialiseCmod(self.nthreads,self.control["cal_source"],wfs.seed)
            #Take correlation image, if don't yet have one, and if doing correlationCentroiding.
            wfs.takeCorrImage(self.control)
            #Take the centroid weighting if we're using it...
            wfs.takeCentWeight(self.control)
            #now calibrate the SHS
            wfs.calibrateSHS(self.control)
            #Now take reference centorids...
            wfs.takeReference(self.control)
            #and reset the output.
            wfs.outputData[:]=0
        if self.fullOutput:
            self.outputData=self.thisObjList[0].wfscentObj.outputData
        else:
            #self.outputData=numpy.zeros((self.thisObjList[0].wfscentObj.nsubaps,2),numpy.float32)
            self.outputData=self.thisObjList[0].wfscentObj.outputData
        #print "wfscent finalInitialisation done"
    def prepareNextIter(self):
        """prepare resource shared object for next computation..."""
        this=self.thisObjList[self.currentIdObjCnt]
        for attr in dir(this):
            if attr not in ["__doc__","__module__","__init__","idstr"]:
                val=getattr(this,attr)
                setattr(self,attr,val)
    def endThisIter(self):
        """Copy things back to the this object..."""
        if self.fullOutput:
            self.outputData=self.wfscentObj.outputData
        else:
            numpy.take(self.wfscentObj.outputData.ravel(),self.wfscentObj.indices,out=self.outputData.ravel())
                
    
    def generateNext(self,msg=None):
        """WFS main iteration loop.
        Not expecting any msgs.
        """
        t1=time.time()
        if self.debug!=None:
            print("INFORMATION:wfscent: Doing generateNext (debug={0:s})".format(
                  str(self.debug)) )
##(old)            print "wfs: Doing generateNext (debug=%s)"%str(self.debug)
        current=self.control["parentType"]
        if self.generate==1:
            if self.newDataWaiting:#this is always 1
                if self.parent[current].dataValid==1:
                    #self.inputData=self.parent.outputData
                    self.inputInvalid=0
                else:
                    print("INFORMATION:wfscent: waiting for data from dm, but "+
                        "not valid")
                    self.inputInvalid=1
                    self.dataValid=0
            if self.inputInvalid==0: # there was an input, so we can integrate...
                wfs=self.wfscentObj # has been copied from thisObjList before generateNext is called...
                if wfs.inputImage is not None:
                    #input data is an image, rather than phase.
                    self.dataValid=1
                    if self.control["zeroOutput"]:
                        wfs.outputData[:]=0
                    else:
                        wfs.runSlopeCalc(self.control)
                        if wfs.subtractTipTilt==-1 or (
                                wfs.subtractTipTilt==1 and self.control["cal_source"]==0 and self.imageOnly==0):
                            N=wfs.nsubaps
                            # subtract average x centroid:
                            wfs.outputData[:,:,0]-=wfs.outputData[:,:,0].sum()/N
                            # subtract average y centroid:
                            wfs.outputData[:,:,1]-=wfs.outputData[:,:,1].sum()/N
                        
                else:
                    self.dataValid=0
                    if wfs.integstepFn!=None and wfs.texp==0 and self.control["useCmod"]:
                        wfs.updateIntegTime(wfs.integstepFn())
                    if wfs.texp<wfs.integtime:
                        # Still integrating...
                        #t=time.time()
                        # Stack up the phases
                        wfs.reorder(self.parent[current].outputData,int(wfs.texp/wfs.tstep))              
                        #print "wfs: Reorder time %g"%(time.time()-t)

                    wfs.texp+=wfs.tstep
                    if wfs.texp>=wfs.integtime+wfs.latency:  # Exposure Complete
                        wfs.texp=0.
                        if self.control["zeroOutput"]:
                            wfs.outputData[:]=0
                        else:
                            wfs.runCalc(self.control)
                            # this should be used for LGS sensors:
                            if wfs.subtractTipTilt==-1 or (
                                wfs.subtractTipTilt==1 and self.control["cal_source"]==0 and self.imageOnly==0):
                                N=wfs.nsubaps
                                # subtract average x centroid:
                                wfs.outputData[:,:,0]-=wfs.outputData[:,:,0].sum()/N
                                # subtract average y centroid:
                                wfs.outputData[:,:,1]-=wfs.outputData[:,:,1].sum()/N
                        self.dataValid=1

                if self.timing:
                    print("INFORMATION:wfscent: time:{0:s}".format( str(time.time()-t1) ))
        else:
            self.dataValid=0
        if self.debug!=None:
            print("INFORMATION:wfscent: Done generateNext (debug={0:s})".format(
                  str(self.debug)) )
        self.generateNextTime=time.time()-t1




    def newCorrRef(self):
        """Grabs current SHS images and sets these as the correlation reference.  Then computes new reference slopes too"""
        cs=self.control["cal_source"]
        self.control["cal_source"]=1
        this=self.thisObjList[0]
        wfs=this.wfscentObj
        wfs.corrPattern=None
        data=wfs.takeCorrImage(self.control)
        #now calibrate the SHS
        #wfs.calibrateSHS(self.control)
        #Now take reference centorids...
        wfs.takeReference(self.control)
        self.control["cal_source"]=cs

        return data




    def drawCents(self,fromcent=0,objNumber=None,mask=None):
        """Draw centroids in a format that can be easily displayed.
        If fromcent is 1, will construct an image from centroid locations,
        otherwise, will use the noisy image (non-fpga only).
        objNumber specifies which resource sharing object to use.  None means
        the one currently in use.
        If mask is set, parts of the images that aren't used in centroid computation will be masked out (ncen).
        If mask is 1, they are masked at max value, and if -1, masked at zero.
        """
        if mask is None:
            mask=self.imgmask
        if objNumber==None:
            try:
                wfsobj=self.wfscentObj
            except:
                print("ERROR:wfscent: getting wfscentObj - assuming "+
                     "thisObjList[0].wfscentObj")
                wfsobj=self.thisObjList[0].wfscentObj
        else:
            wfsobj=self.thisObjList[objNumber].wfscentObj
        
        if fromcent:
            pow=0.125
            divisor=2*((wfsobj.nimg-1)**2)**pow
            self.shimg[:,]=0.
            for i in xrange(wfsobj.nsubx):
                for j in xrange(wfsobj.nsubx):
                    xx=wfsobj.centx[i,j]+wfsobj.nimg/2-0.5
                    yy=wfsobj.centy[i,j]+wfsobj.nimg/2-0.5
                    xf=int(xx)
                    yf=int(yy)
                    xc=int(math.ceil(xx))
                    yc=int(math.ceil(yy))
                    self.shimg[xf+i*wfsobj.nimg,yf+j*wfsobj.nimg]=(1-(xx-xf))*(1-(yy-yf))
                    if xf!=xc:
                        self.shimg[xc+i*wfsobj.nimg,yf+j*wfsobj.nimg]=(1-(xx-xc))*(1-(yy-yf))
                        if yf!=yc:
                            self.shimg[xc+i*wfsobj.nimg,yc+j*wfsobj.nimg]=(1-(xx-xc))*(1-(yy-yc))

                    if yf!=yc:
                        self.shimg[xf+i*wfsobj.nimg,yc+j*wfsobj.nimg]=(1-(xx-xf))*(1-(yy-yc))
        else:
            if self.control["useCmod"]:
                bimg=wfsobj.cmodbimg
            else:
                bimg=wfsobj.bimg.astype("f")
            for i in xrange(wfsobj.nsubx):     # Loop over subaps
                for j in xrange(wfsobj.nsubx):
                    self.shimg[i*wfsobj.nimg:(i+1)*wfsobj.nimg,j*wfsobj.nimg:(j+1)*wfsobj.nimg]=\
                        (bimg[i,j]*wfsobj.subflag[i,j])#.astype("f")    # Tessalate up for WFS display
        result=self.shimg[:wfsobj.nimg*wfsobj.nsubx,:wfsobj.nimg*wfsobj.nsubx]
        if mask!=0:
            if mask==1:
                maskval=max(result.flat)
            else:
                maskval=0
            ncen=wfsobj.ncen
            nimg=wfsobj.nimg
            if ncen<nimg:
                s=(nimg-ncen)/2#starting point
                e=nimg-s#end point
                for i in xrange(wfsobj.nsubx):
                    result[i*nimg:i*nimg+s]=maskval
                    result[i*nimg+e:i*nimg+nimg]=maskval
                    result[:,i*nimg:i*nimg+s]=maskval
                    result[:,i*nimg+e:i*nimg+nimg]=maskval

        return result
    
    def drawCorrelation(self,img,objno=None,mask=None):
        """img can be wfsobj.corrimg, or wfsobj.corrPatternUser or "corrimg" or "corrPattern"
        img.shape==nsubx,nsubx,nimg,nimg
        """
        if mask is None:
            mask=self.imgmask
        if objno==None:
            wfsobj=self.wfscentObj
        else:
            wfsobj=self.thisObjList[objno].wfscentObj
        if img=="corrimg":
            img=wfsobj.corrimg
        elif img=="corrPattern":
            img=wfsobj.corrPatternUser
        nimg=img.shape[2]
        for i in xrange(wfsobj.nsubx):
            for j in xrange(wfsobj.nsubx):
                 # Tessalate up for WFS display:
                self.shimg[i*nimg:(i+1)*nimg,j*nimg:(j+1)*nimg]=img[i,j]
        result=self.shimg[:nimg*wfsobj.nsubx,:nimg*wfsobj.nsubx]
        if mask!=0:
            if mask==1:
                maskval=max(result.flat)
            else:
                maskval=0
            ncen=wfsobj.ncen
            #nimg=wfsobj.nimg
            if ncen<nimg:
                s=(nimg-ncen)/2#starting point
                e=nimg-s#end point
                for i in xrange(wfsobj.nsubx):
                    result[i*nimg:i*nimg+s]=maskval
                    result[i*nimg+e:i*nimg+nimg]=maskval
                    result[:,i*nimg:i*nimg+s]=maskval
                    result[:,i*nimg+e:i*nimg+nimg]=maskval
        return result
            

    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        txt=""
        this=self.thisObjList[self.sentPlotsCnt]
        if this.idstr==None or this.idstr=="":
            id=""
        else:
            id=" (%s)"%this.idstr
        if self.sentPlotsCnt==0:
            #outputData is only valid for one object at a time, when that has just run...
            txt+="""<plot title="WFS SH img%s (the SHS images)" cmd="data=%s.drawCents(0)" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            txt+="""<plot title="Change mask (-1,0,1) %s (to change the colour of guard pixels)" cmd="data=%s.imgmask=(%s.imgmask+2)%%3-1" when="cmd" ret="data" texttype="1"/>"""%(id,objname,objname)
            if self.imageOnly==0:
                txt+="""<plot title="XCentroids%s (display X slopes)" cmd="data=%s.wfscentObj.outputData[:,:,0]" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
                txt+="""<plot title="YCentroids%s (display Y slopes)" cmd="data=%s.wfscentObj.outputData[:,:,1]" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
                txt+="""<plot title="1D centroids%s (display slopes as a 1D plot)" cmd="data=%s.wfscentObj.outputData.ravel()" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            elif self.imageOnly==1:
                txt+="""<plot title="Centroids%s (display slopes as a 1D plot)" cmd="data=%s.wfscentObj.outputData.ravel()" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            else:
                txt+="""<plot title="Centroids%s (display slopes as a 1D plot)" cmd="data=%s.wfscentObj.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            txt+="""<plot title="Centroids 2D%s (a 2D representation of slopes)" cmd="data=%s.drawCents(1)" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            for i in range(len(self.thisObjList)):
                this=self.thisObjList[i]
                wfs=this.wfscentObj
                if type(this.laserGuideStar)!=type(None):
                    if type(this.laserGuideStar)==numpy.ndarray:
                        #txt+="""<plot title="LGS elongation%s" cmd="data=%s.thisObjList[%d].laserGuideStar" ret="data" type="pylab" when="cmd" palette="gray"/>\n"""%(id,objname,i)
                        txt+="""<plot title="LGS elongation%s (show the LGS PSF)" cmd="data=%s.thisObjList[%d].wfscentObj.reformatImg(%s.thisObjList[%d].laserGuideStar)" ret="data" type="pylab" when="cmd" palette="gray"/>\n"""%(id,objname,i,objname,i)
                    else:
                        txt+="""<plot title="LGS elongation%s (show the PSF PSF)" cmd="data=%s.thisObjList[%d].laserGuideStar.subapImage" ret="data" type="pylab" when="cmd" palette="gray"/>\n"""%(id,objname,i)
            for i in range(len(self.thisObjList)):
                this=self.thisObjList[i]
                wfs=this.wfscentObj
                if wfs.correlationCentroiding:
                    txt+="""<plot title="Correlation%s (show the current correlation image)" cmd="data=%s.drawCorrelation('corrimg')" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
                    txt+="""<plot title="Correlation PSF%s (show the correlation reference image)" cmd="data=%s.drawCorrelation('corrPattern',mask=0)" ret="data" type="pylab" when="cmd" palette="gray"/>"""%(id,objname)
                    txt+="""<plot title="Get new corr ref%s" cmd="data=%s.newCorrRef()" ret="data" type="pylab" when="cmd" palette="gray"/>"""%(id,objname)

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
        paramList.append(base.dataType.dataType(description="dataType",typ="eval",val="this.globals.fpDataType",comment="Array numpy data type"))
        paramList.append(base.dataType.dataType(description="timing",typ="i",val="0",comment="Timing information"))
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="ntel",typ="eval",val="this.globals.npup",comment="Pixels for telescope"))
        paramList.append(base.dataType.dataType(description="telSec",typ="f",val="8.",comment="TODO: Telescope secondary diameter (m)"))
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam,wfs_nsubx,wfs_minarea)",comment="Telescope pupil"))
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))
        paramList.append(base.dataType.dataType(description="tstep",typ="f",val="0.005",comment="TODO: timestep."))        
        return paramList


