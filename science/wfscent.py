"""A module for simulating a WFS and centroider, this can make use of the FPGAs, offering speedups of up to 400 times."""
#$Id: wfscent.py,v 1.65 2011/12/06 05:14:39 ali Exp $
### Shack Hartmann WFS simulation object  #####
# Integrate WFS CCD images

import math
import numpy,numpy.random
import cmod.imgnoise
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
        self.imageOnly=self.config.getVal("imageOnly",default=0)# 0 to return slopes, 
                                                                # 1 to return image as nsubx,nsubx,nimg,nimg, 
                                                                # 2 to return image as a 2d image.
        if forGUISetup==1:#this won't be set if resource sharing - whats the point!  So, this is ok.
            nsubx=self.config.getVal("wfs_nsubx")
            nimg=self.config.getVal("wfs_nimg")
            fullOutput=self.config.getVal("fullWFSOutput",default=1)
            if fullOutput==1:
                if self.imageOnly==0:
                    self.outputData=[(nsubx,nsubx,2),numpy.float32]
                elif self.imageOnly==1:
                    self.outputData=[(nsubx,nsubx,nimg,nimg),numpy.float32]
                else:
                    self.outputData=[(nsubx*nimg,nsubx*nimg),numpy.float32]

            else:
                n=self.config.getVal("wfs_n")                                  
                pupfn=self.config.getVal("pupil")
                if type(pupil)!=numpy.ndarray:
                    pupfn=pupfn.fn
                pupfn=pupfn.astype(numpy.float32)
                wfs_minarea=self.config.getVal("wfs_minarea")
                frac=wfs_minarea*n*n
                nsubaps=0
                for i in xrange(nsubx):        
                    for j in xrange(nsubx):
                        if pupfn[i*n:(i+1)*n,j*n:(j+1)*n].sum()>frac:
                            nsubaps+=1
                self.outputData=[(nsubaps,2),numpy.float32]
        else:
            self.debug=debug
            self.fpid=None
            self.fullOutput=self.config.getVal("fullWFSOutput",default=1)
            self.fpDataType=self.config.getVal("fpDataType",default=numpy.float32)
            self.doneFinalInit=0
##             self.atmosPhaseType=self.config.getVal("atmosPhaseType",default="phaseonly")
##             if self.atmosPhaseType not in ["phaseonly","phaseamp","realimag"]:
##                 raise Exception("wfscent: atmosPhaseType not known %s"%self.atmosPhaseType)
            #seed here - no point seeding each resource sharing object...
            numpy.random.seed(config.getVal("wfs_seed",default=None,raiseerror=0))
            #random.seed(config.getVal("wfs_seed",default=0))
            cmod.imgnoise.seed(config.getVal("wfs_rand_seed",default=0))
            if args.has_key("fpprecision"):
                if args["fpprecision"]==32:
                    self.fpDataType=numpy.float32
                else:
                    self.fpDataType=numpy.float64
            self.sentPlotsCnt=0
#             self.initFPGA=self.config.getVal("initFPGA",default=0)#whether to initialise the FPGA... (load binary, set up arrays etc).
#             self.FPGABitFile=None
#             self.ignoreFPGAOpenFailure=0
#             if self.initFPGA:
#                 self.FPGABitFile=self.config.getVal("FPGAWFSBitFile",default=string.join(__file__.split("/")[:-2]+["fpga","wfspipe.bin.ufp"],"/"))
#                 self.ignoreFPGAOpenFailure=self.config.getVal("ignoreFPGAOpenFailure",default=0)
# ##             self.waitFPGA=self.config.getVal("waitFPGA",1)
#                 useFPGA=self.config.getVal("useFPGA",default=0)#whether to use the FPGA initially (global)
#             else:
#                 useFPGA=0
#            useCell=self.config.getVal("useCell",default=0)
            useCmod=self.config.getVal("useCmod",default=1)
            self.printLinearisationForcing=self.config.getVal("printLinearisationForcing",default=0)
            #self.cellseed=self.config.getVal("cellseed",default=1)
            self.cmodcentseed=self.config.getVal("cmodcentseed",default=0)
            self.nthreads=self.config.getVal("nthreads",default="all")#usually an integer... or "all"
            if self.nthreads=="all":#use all available CPUs...
                self.nthreads=self.config.getVal("ncpu")#getCpus()
                print "wfscent: Using %d threads"%self.nthreads
            calSource=self.config.getVal("calsource")
            self.imgmask=1#value used when creating sh img in drawCents()
#            self.control={"cal_source":calSource,"useFPGA":useFPGA,"useCell":useCell,"useCmod":useCmod,"zeroOutput":0,"parentType":"closed"}#ZeroOutput can be used when calibrating a pseudo-open loop system.
            self.control={"cal_source":calSource,"useCmod":useCmod,"zeroOutput":0,"parentType":"closed"}#ZeroOutput can be used when calibrating a pseudo-open loop system.
            self.lastCalSource=0#resourcesharing - todo...?

##             # Telescope pupil phase array size
##             self.npup=self.config.getVal("npup")                                    
##             # No. of subaps across tel. pupil
##             self.wfs_nsubx=self.config.getVal("wfs_nsubx")                          
##             # Phase array size for subaps (phasesize)
##             self.wfs_n=self.config.getVal("wfs_n")                                  
##             # FFT size for subap image calc (fftsize)
##             self.wfs_nfft=self.config.getVal("wfs_nfft")                            
##             # Image size for subap (pixels) (ie size after binning).
##             self.wfs_nimg=self.config.getVal("wfs_nimg")                            
##             self.bin=0
##             if(self.wfs_nfft!=self.wfs_nimg):
##                 self.bin=1
##                 if Numeric.floor(self.wfs_nfft/self.wfs_nimg)!=Numeric.ceil(self.wfs_nfft/self.wfs_nimg):
##                     print "WARNING: Non-integer binning of image - software (binimgmodule.so) cannot handle this (though FPGA will be fine!)."


##             self.xtiltfn=((Numeric.fromfunction(self.tilt,(self.wfs_n,self.wfs_n))-float(self.wfs_n)/2.+0.5)/float(self.wfs_n)).astype(self.fpDataType)# subap tilt fn
##             self.ytiltfn=Numeric.transpose(self.xtiltfn)

##             self.subimg=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.wfs_nfft,self.wfs_nfft),Numeric.Float64)      # SH sub-images (high LL)
##             self.bimg=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.wfs_nimg,self.wfs_nimg),Numeric.Float64)
            #self.outputData=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.wfs_nimg,self.wfs_nimg),Numeric.Int32)#same dimensions as bimg, different type.

##             if self.fpDataType==Numeric.Float64:
##                 self.fftwPlan=cmod.mkimg.setup(self.subimg[0][0])                   # Setup imaging FFTs for subaps
##             else:
##                 self.fftwPlan=cmod.mkimgfloat.setup(self.subimg[0][0])

##             self.texp=0.                                                            # Elapsed exposure time
##             self.tstep=self.config.getVal("tstep")                                  # Simulation tick length
##             self.wfs_int=self.config.getVal("wfs_int")                              # WFS integration time
##             self.wfs_read_mean=self.config.getVal("wfs_read_mean")                  # WFS Readnoise e-
##             self.wfs_read_sigma=self.config.getVal("wfs_read_sigma")
##             self.fpga_read_sigma=self.wfs_read_sigma/26.11
##             self.fpga_read_mean=self.wfs_read_mean-self.fpga_read_sigma*127
            self.timing=self.config.getVal("timing")
##             self.wfs_lat=self.config.getVal("wfs_lat")                              # WFS readout latency
##             self.skybrightness=self.config.getVal("wfs_skybrightness")
##             self.shimg=Numeric.zeros((self.wfs_nsubx*self.wfs_nimg,self.wfs_nsubx*self.wfs_nimg),Numeric.Float64)# Tessalated SH image for display
            #self.cal_source=config.getVal("calsource")
##             self.tmpImg=Numeric.zeros((self.wfs_nfft,self.wfs_nfft),self.fpDataType)
            #number of integrations before wfs is read out.
##             self.nIntegrations=int(math.ceil(self.wfs_int/self.tstep))
            #time taken to DMA memory to FPGA (estimate)
##             self.waitFPGATime=self.config.getVal("waitFPGATime",default=(self.wfs_nsubx*self.wfs_n)**2*self.nIntegrations*5e-9)
##             if math.ceil(self.wfs_int/self.tstep)!=math.floor(self.wfs_int/self.tstep):
##                 print "Warning: wfs - Integration times is not a whole number of timesteps - you might misinterpret the results... %g %g"%(self.wfs_int,self.tstep)

##             self.subflag=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx),Numeric.Int8)
##             if self.subflag.itemsize()==8:
##                 print "WARNING: untested with 8 byte longs...(wfs)"
##             self.subarea=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx),Numeric.Float64)
##             self.pupsub=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.wfs_n,self.wfs_n),self.fpDataType)
##             self.pupil=self.config.getVal("pupil")
##             wfs_minarea=self.config.getVal("wfs_minarea")#0.5...                                      # Min unvignetted subap area to use - why here ????
##             n=self.wfs_n
##             for i in range(self.wfs_nsubx):        
##                 for j in range(self.wfs_nsubx):
##                     self.pupsub[i][j]=self.pupil.fn[i*n:(i+1)*n,j*n:(j+1)*n].astype(self.fpDataType)    # Get pupil fn over subaps
##                     self.subarea[i][j]=Numeric.sum(Numeric.sum(self.pupsub[i][j]))                          
##                     if(self.subarea[i][j]>(wfs_minarea*n*n)):# Flag vignetted subaps 
##                         self.subflag[i][j]=1
                        
##             self.dmpupil=Numeric.zeros((self.npup,self.npup),Numeric.Float64)                 # DM pupil-set to zero for 
##             for i in range(self.wfs_nsubx):                                                  # flagged subaps: assume they are    
##                 for j in range(self.wfs_nsubx):                                              # tilted out of science fov
##                     if(self.subflag[i][j]==1):
##                         self.dmpupil[i*n:(i+1)*n,j*n:(j+1)*n]=1.





##             self.wfs_bandwidth=self.config.getVal("wfs_bandwidth")                  # WFS Optical bandwidth
##             self.wfs_thruput=self.config.getVal("wfs_thruput")                      # Thruput to the WFS
##             self.wfs_phot_rate_factor=self.config.getVal("wfs_phot_rate_factor")    # Photons/sec/area/angstrom
##             #self.wfs_int=self.config.getVal("wfs_int")                              
##             self.wfs_mag=self.config.getVal("wfs_mag")                              # Guide star magnitude
##             self.telDiam=self.config.getVal("telDiam")                              # Telescope aperture diameter
##             self.sig=self.wfs_sig()                                                 # calculate WFS signal (e-/subap/integn)
            #The FPGA version of this is in 20.7 format, and also should be reduced by the number of phase pixels per subap...
##             self.sigFPGA=int(self.sig*2**7/(self.wfs_n**2))&0x7ffffff#20.7 format.
##             if int(self.sig*2**7)>0x7ffffff:#dont divide by wfs_n^2 here.
##                 print "wfscent: Warning - SIG too large, will overflow in FPGA"

##             self.wfs_ncen=self.config.getVal("wfs_ncen")            # Centroiding box size (pixels)
            #if self.wfs_ncen!=self.wfs_nimg:
            #    print "WARNING: ncen!=nimg, FPGA will ignore ncen value and use whole box for centroiding."
##             self.wfs_floor=self.config.getVal("wfs_floor")          # Centroiding floor value
##             wfs_minarea=self.config.getVal("wfs_minarea")           # Min unvignetted subap area to use

##             self.cenmask=Numeric.zeros((self.wfs_nimg,self.wfs_nimg),Numeric.Float32)             # Define centroiding mask
##             self.cenmask[self.wfs_nimg/2-self.wfs_ncen/2:self.wfs_nimg/2+self.wfs_ncen/2,self.wfs_nimg/2-self.wfs_ncen/2:self.wfs_nimg/2+self.wfs_ncen/2]=1.
##             self.tilt_indx = (Numeric.array(range(self.wfs_nimg),Numeric.Float64))-float(self.wfs_nimg/2)+0.5#Index fns for centroiding
##             self.centDisplay=Numeric.zeros((self.wfs_nsubx*self.wfs_nimg,self.wfs_nsubx*self.wfs_nimg),Numeric.Float32)
##             self.fpgaInitialised=0
            
##             if self.initFPGA and self.testFPGAUsage():
##                 self.initialiseFPGA()
##             if self.fpgaInitialised==0:
##                 self.initDummyFPGA()#set up some memory...
##                 self.canUseFPGA=0
            for i in xrange(len(self.idstr)):
                idstr=self.idstr[i]
                parent=self.parentList[i]
                self.initialise(parent,idstr)

##     def getCpus(self):
##         ncpu=None
##         try:
##             lines=open("/proc/cpuinfo").readlines()
##             ncpu=0
##             for line in lines:
##                 if ("processor" in line) and (":" in line):
##                     ncpu+=1
##             if ncpu==0:
##                 print "Warning - couldn't determine number of CPUs, assuming 1"
##                 ncpu=1
##         except:
##             pass
##         if ncpu==None:
##             try:
##                 cmd="sysctl -a hw | grep ncpu | tail -n 1"
##                 print "/proc/cpuinfo not found - may be OSX?  Trying commandline '%s'"%cmd
##                 line=popen2.popen2(cmd)[0].read().strip()
##                 ncpu=int(line.split("=")[1])
##             except:
##                 ncpu=1
##                 print "Cannot detemine number of CPUs - assuming 1"
##         return ncpu
    def initialise(self,parent,idstr):
        """Perform initialisation for this resource-sharing object."""
        if type(parent)!=type({}):
            parent={"closed":parent}
        this=base.aobase.resourceSharer(parent,self.config,idstr,self.moduleName)
        print "Initialising centroid object %s %s"%(str(parent),str(idstr))
        self.thisObjList.append(this)
        #fpDataType=this.config.getVal("fpDataType",default=Numeric.Float32) # Numeric.Float64
        #if self.fpDataType!=None:#specified by args{}
        #    fpDataType=self.fpDataType
        atmosPhaseType=this.config.getVal("atmosPhaseType",default="phaseonly")
        if atmosPhaseType not in ["phaseonly","phaseamp","realimag"]:
            raise Exception("wfscent: atmosPhaseType not known %s"%atmosPhaseType)
        # useFPGA=self.config.getVal("useFPGA",default=0)#whether to use the FPGA (shouldn't change)
        # if useFPGA:
        #     waitFPGA=self.config.getVal("waitFPGA",default=1)
        # else:
        #     waitFPGA=1
        # Telescope pupil phase array size
        #npup=this.config.getVal("npup")                                    
        # No. of subaps across tel. pupil
        wfs_nsubx=this.config.getVal("wfs_nsubx")                          
        # Phase array size for subaps (phasesize)
        wfs_n=this.config.getVal("wfs_n")                                  
        # FFT size for subap image calc (fftsize)
        wfs_nfft=this.config.getVal("wfs_nfft",default=wfs_n*2)
        clipsize=this.config.getVal("wfs_clipsize",default=wfs_nfft)
        # Image size for subap (pixels) (ie size after binning).
        wfs_nimg=this.config.getVal("wfs_nimg",default=clipsize/2)
        #todo: this bit should really be in the centroid module. Not needed here...
        bin=0
        if(clipsize!=wfs_nimg):
            bin=1
            if (float(clipsize)/wfs_nimg)%1!=0:
                print "WARNING: wfscent - Non-integer binning of image - software (binimgmodule.so) cannot handle this (though FPGA will be fine!)."
        tstep=this.config.getVal("tstep")                                  # Simulation tick length
        wfs_int=this.config.getVal("wfs_int")                              # WFS integration time
        wfs_read_mean=this.config.getVal("wfs_read_mean")                  # WFS Readnoise e-

        wfs_read_sigma=this.config.getVal("wfs_read_sigma")
        threshType=this.config.getVal("threshType",default=0)
        #todo: put in centroid module (next 2 lines).
        #self.fpga_read_sigma=self.wfs_read_sigma/26.11
        #self.fpga_read_mean=self.wfs_read_mean-self.fpga_read_sigma*127
        wfs_lat=this.config.getVal("wfs_lat")                              # WFS readout latency
        skybrightness=this.config.getVal("wfs_skybrightness")
        nIntegrations=int(math.ceil(wfs_int/tstep))        
        # if useFPGA:
        #     waitFPGATime=this.config.getVal("waitFPGATime",default=(wfs_nsubx*wfs_n)**2*nIntegrations*5e-9)
        # else:
        #     waitFPGATime=0
        if (wfs_int/tstep)%1!=0:
            print "Warning: wfs - Integration times is not a whole number of timesteps - you might misinterpret the results... %g %g"%(wfs_int,tstep)
        pupil=this.config.getVal("pupil")
        wfs_minarea=this.config.getVal("wfs_minarea")#0.5... # Min unvignetted subap area to use - why here ????


        
        wfs_bandwidth=this.config.getVal("wfs_bandwidth")                  # WFS Optical bandwidth
        wfs_thruput=this.config.getVal("wfs_thruput")                      # Thruput to the WFS
        wfs_phot_rate_factor=this.config.getVal("wfs_phot_rate_factor")    # Photons/sec/area/angstrom
        wfs_mag=this.config.getVal("wfs_mag")                              # Guide star magnitude
        telDiam=this.config.getVal("telDiam")                              # Telescope aperture diameter
        wfs_ncen=this.config.getVal("wfs_ncen",default=wfs_nimg)            # Centroiding box size (pixels)
        wfs_floor=this.config.getVal("wfs_floor")          # Centroiding floor value
        sig=this.config.getVal("wfs_sig",default=None,raiseerror=0)
        if sig==None:
            sig=util.centroid.wfs_sig(wfs_bandwidth,wfs_thruput,wfs_phot_rate_factor,telDiam,wfs_nsubx,wfs_int,wfs_mag)
        this.laserGuideStar=None
        addPoisson=this.config.getVal("addPoisson",default=1)
        #Now look to see if we're using LGS...
        laserGuideStar=this.config.getVal("laserGuideStar",default=None,raiseerror=0)#this can be an instance of util.guideStar.wfs
        if laserGuideStar==None:
            lgsType=this.config.getVal("lgsType",default=None,raiseerror=0)#can be None, "sodium" or "rayleigh".
            if lgsType in ["sodium","rayleigh"]:
                # some extra configuration is needed for setting up the spot elongation...
                lgsPower=this.config.getVal("lgsPower",default=5.)
                telDiam=this.config.getVal("telDiam")
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
        if laserGuideStar!=None:
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
                print "Generating LGS spot patterns... (if these are all the same and take a while, you might be better specifying this in the param file...)"
                laserGuideStar.wfs_image(off=0.5)
                laserGuideStar.sumSubaps()
                print "Finished generating LGS spot patterns."
                spotpsf=laserGuideStar.spotPsf
                sig=laserGuideStar.sig
                this.laserGuideStar=laserGuideStar
            
        subtractTipTilt=this.config.getVal("subtractTipTilt",default=int(this.laserGuideStar!=None))
        if type(sig)==type(0.):
            print "wfs sig is %g phot/subap"%sig
        else:
            print "wfs sig is array with max %g phot/subap"%max(sig.flat)
        if type(spotpsf)==type(None):
            spotpsf=this.config.getVal("spotpsf",default=None,raiseerror=0)#a spot PSF, eg an airy disc (eg createAiryDisc(self.fftsize,self.fftsize/2,0.5,0.5)), or LGS elongated spots.
        opticalBinning=this.config.getVal("opticalBinning",default=0)#whether to simulate optical binning (beam splitter and 2 cylindrical lenslet arrays and 2 detectors)
        magicCentroiding=this.config.getVal("magicCentroiding",default=0)
        linearSteps=this.config.getVal("linearSteps",default=None,raiseerror=0)
        calNCoeff=this.config.getVal("calNCoeff",default=0)
        stepRangeFrac=this.config.getVal("stepRangeFrac",default=1.)
        phslam=None
        wfslam=None
        atmosGeom=this.config.getVal("atmosGeom",raiseerror=0)
        if atmosGeom!=None:
            wfslam=atmosGeom.sourceLambda(idstr)
            phslam=atmosGeom.phaseLambda(idstr)
        if wfslam==None:
            wfslam=this.config.getVal("sourceLam")
            print "Warning - wfslam not in atmosGeom, using %g"%wfslam
        if phslam==None:
            phslam=wfslam
        phaseMultiplier=phslam/wfslam
        if phaseMultiplier!=1:
            print "Warning - wfscent %s has phaseMultiplier of %g"%(str(idstr),phaseMultiplier)
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
#        this.wfscentObj=util.centroid.centroid(wfs_nsubx,pup=pupil,oversamplefactor=None,readnoise=wfs_read_sigma,readbg=wfs_read_mean,addPoisson=1,noiseFloor=wfs_floor,binfactor=None,sig=sig,skybrightness=skybrightness,warnOverflow=None,atmosPhaseType=atmosPhaseType,fpDataType=self.fpDataType,useFPGA=useFPGA,waitFPGA=waitFPGA,waitFPGATime=waitFPGATime,phasesize=wfs_n,fftsize=wfs_nfft,clipsize=clipsize,nimg=wfs_nimg,ncen=wfs_ncen,tstep=tstep,integtime=wfs_int,latency=wfs_lat,wfs_minarea=wfs_minarea,spotpsf=spotpsf,opticalBinning=opticalBinning,useCell=self.control["useCell"],waitCell=1,usecmod=self.control["useCmod"],subtractTipTilt=subtractTipTilt,magicCentroiding=magicCentroiding,linearSteps=linearSteps,stepRangeFrac=stepRangeFrac,phaseMultiplier=phaseMultiplier,centWeight=centWeight,correlationCentroiding=correlationCentroiding,corrThresh=corrThresh,corrPattern=corrPattern,threshType=threshType,imageOnly=self.imageOnly,calNCoeff=calNCoeff,useBrightest=useBrightest)
        this.wfscentObj=util.centroid.centroid(wfs_nsubx,pup=pupil,oversamplefactor=None,readnoise=wfs_read_sigma,readbg=wfs_read_mean,addPoisson=addPoisson,noiseFloor=wfs_floor,binfactor=None,sig=sig,skybrightness=skybrightness,warnOverflow=None,atmosPhaseType=atmosPhaseType,fpDataType=self.fpDataType,phasesize=wfs_n,fftsize=wfs_nfft,clipsize=clipsize,nimg=wfs_nimg,ncen=wfs_ncen,tstep=tstep,integtime=wfs_int,latency=wfs_lat,wfs_minarea=wfs_minarea,spotpsf=spotpsf,opticalBinning=opticalBinning,usecmod=self.control["useCmod"],subtractTipTilt=subtractTipTilt,magicCentroiding=magicCentroiding,linearSteps=linearSteps,stepRangeFrac=stepRangeFrac,phaseMultiplier=phaseMultiplier,centWeight=centWeight,correlationCentroiding=correlationCentroiding,corrThresh=corrThresh,corrPattern=corrPattern,threshType=threshType,imageOnly=self.imageOnly,calNCoeff=calNCoeff,useBrightest=useBrightest,printLinearisationForcing=self.printLinearisationForcing)


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
#        fpgarequired=0
        nsubxmax=0
        nsubxnfftmax=0
        maxarrsize=0
        nIntMax=0
        maxnpup=0
        rpmax=0
        imgsizemax=0
        atmosfactor=1#phaseonly...
#        fpgaObj=None
#        fpgaarr=None
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
            if wfs.nsubx*wfs.fftsize>nsubxnfftmax:
                nsubxnfftmax=wfs.nsubx*wfs.fftsize
            if wfs.nsubx*wfs.phasesize>maxnpup:
                maxnpup=wfs.nsubx*wfs.phasesize
            if wfs.nsubx*wfs.nsubx*wfs.nIntegrations*wfs.phasesize*wfs.phasesize>rpmax:
                atmosPhaseType=this.config.getVal("atmosPhaseType",default="phaseonly")
                print this.config.searchOrder
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
#            if wfs.canUseFPGA==1:
#                fpgarequired=1
#                fpgaObj=wfs
            arrsize=wfs.nsubx*wfs.nsubx*wfs.nIntegrations*wfs.phasesize*wfs.phasesize*4*atmosfactor+wfs.nsubx*wfs.nsubx*2*4
            if arrsize>maxarrsize:
                maxarrsize=arrsize
        #now allocate memories...
        reorderedPhsMem=None
        # wfs=fpgaObj#any object that can use the FPGAs
        # if fpgarequired:#use the first object to do some generic setup...
        #     fpid,fpgaInfo=wfs.initialiseFPGA(ignoreFailure=self.ignoreFPGAOpenFailure,fpgaBitFile=self.FPGABitFile)#load the fpga binary.
        #     fpgaarr=wfs.setupFPGAArray(fpid,maxarrsize)#allocate FPGA buffer
        # else:
#        fpid=None
#        fpgaInfo=None
        if canSharePhs:
            reorderedPhsMem=numpy.zeros((rpmax,),numpy.float32)
        if self.imageOnly:
            outputDataMem=numpy.zeros((imgsizemax*imgsizemax,),numpy.float32)
        else:
            outputDataMem=numpy.zeros((nsubxmax,nsubxmax,2),numpy.float32)
        bimgMem=numpy.zeros((imgsizemax*imgsizemax,),numpy.float64)#needs 64 bits for cmod.binimg...
        self.shimg=numpy.zeros((imgsizemax,imgsizemax),numpy.float32)
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
            #share the fpga ID and info.
            #wfs.initialiseFPGA(fpid=fpid,fpgaInfo=fpgaInfo,ignoreFailure=self.ignoreFPGAOpenFailure)

            #do we really want to initCell for every wfscentObj?  No.
            #wfs.initialiseCell(nspu=6,calsource=0,showCCDImg=0,allCents=1,cellseed=1)
            wfs.finishInit()
            #wfs.initialiseCell(nspu=6,calsource=self.control["cal_source"],showCCDImg=0,allCents=1,cellseed=self.cellseed)
            wfs.initialiseCmod(self.nthreads,self.control["cal_source"],self.cmodcentseed)
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
                
##     def testFPGAUsage(self):
##         """Checks variables are suitable for FPGA use"""
##         self.canUseFPGA=haveFPGA
##         if self.atmosPhaseType!="phaseonly":
##             print "WARNING: Cannot use FPGA - atmosPhaseType must be phaseonly"
##             self.canUseFPGA=0
##         if self.fpDataType!=Numeric.Float32:
##             print "WARNING: Cannot use FPGA - fpDataType must be Float32"
##             self.canUseFPGA=0
##         if self.wfs_nfft not in [8,16,32]:
##             print "WARNING: FPGA cannot use this FFT array size"
##             self.canUseFPGA=0
##         if self.nIntegrations<1 or self.nIntegrations>63:
##             print "WARNING: Illegal number of integrations for FPGA - must be less than 64"
##             self.canUseFPGA=0
##         if self.wfs_nimg<2 or self.wfs_nimg>32:
##             print "WARNING: Illegal pixel size for centroiding in FPGA (wfs_nimg) - must be 2-32"
##             self.canUseFPGA=0
##         if self.wfs_n<1 or self.wfs_n>32:
##             print "WARNING: Illegal phase size for use in FPGA (wfs_n) - must be <32"
##             self.canUseFPGA=0
##         if int(self.sig*2**7)>0x7ffffff:
##             print "WARNING: Signal is too bright for use in FPGA - should be less than 0xfffff"
##             self.canUseFPGA=0
##         if self.skybrightness>0xffff:
##             print "WARNING: Sky background is too bright for use in FPGA - should be less than 0xffff"
##             self.canUseFPGA=0
##         if self.wfs_floor>0xffff:
##             print "WARNING: WFS floor (threshold) value is too high for use in FPGA - should be less than 0xffff"
##             self.canUseFPGA=0
##         if self.fpga_read_mean>0xffff or self.fpga_read_sigma>0xff:
##             print "WARNING: CCD readout noise is too high for use in FPGA (mean or sigma)"
##             self.canUseFPGA=0
##         if self.wfs_nsubx<1 or self.wfs_nsubx>1023:
##             print "WARNING: Number of subapertures is too large for use in FPGA (needs nsubx<1024)"
##             self.caUseFPGA=0
##         return self.canUseFPGA

##     def initialiseFPGA(self):
##         """Set up the FPGA memory etc, load the FPGA bitstream etc..."""
##         if self.fpid==None:
##             try:                
##                 self.fpid=fpga.open(reset=0,start=0)
##             except:
##                 if self.ignoreFPGAOpenFailure:
##                     self.fpgaInitialised=0
##                     print "WARNING: wfscent - failed to initialise FPGA"
##                     return
##                 else:
##                     raise
##         fpga.load(self.fpid,self.FPGABitFile)
##         fpga.reset(self.fpid)
##         time.sleep(0.001)
##         fpga.start(self.fpid)
##         time.sleep(0.001)
##         fpga.writeReg(self.fpid,0x2,6)#stop fpga pipeline

        # now set up the FPGA DMA arrays...
##         if self.atmosPhaseType=="phaseonly":
##             atmosfactor=1
##         else:
##             raise Exception("atmosphasetype...")
##         arrsize=self.wfs_nsubx*self.wfs_nsubx*self.nIntegrations*self.wfs_n*self.wfs_n*4*atmosfactor+self.wfs_nsubx*self.wfs_nsubx*2*4
##         if arrsize<4*1024*1024:
##             arrsize=4*1024*1024#needs 4MB array for loading QDR memory.
##         if arrsize>1024*1024*1024:
##             print "Warning: Pupil pixel size is too large for single FPGA array - will use multiply arrays, but speed will be reduced (FPGA can access only a 1GB buffer)."
##             self.doPartialFPGA=1
##             #now see how many subaps can fit at once...
##             self.fittedSubaps=1024*1024*1024/(self.nIntegrations*self.wfs_n*self.wfs_n*4*atmosfactor+2*4)
##             arrsize=self.fittedSubaps*(self.nIntegrations*self.wfs_n*self.wfs_n*4*atmosfactor+2*4)
##             self.fpgaarr=fpga.mallocHostMem(self.fpid,(arrsize,),Numeric.Int8)
##             #temporary input and output arrays for FPGA to access.
##             if atmosfactor==1:
##                 self.fpgaInArr=cmod.utils.arrayFromArray(self.fpgaarr,(self.fittedSubaps,self.nIntegrations,self.wfs_n,self.wfs_n),Numeric.Float32)
##             else:
##                 pass #never get here cos error has been raised.
            
##             #create the input and output to copy from and to...
##             if self.atmosPhaseType=="phaseonly":
##                 self.reorderedPhs=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n),Numeric.Float32)
##             elif self.atmosPhaseType=="phaseamp":#shouldn't get here, cos error...
##                 self.reorderedPhs=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n,2),Numeric.Float32)
##             else:
##                 raise Exception("atmosPhaseType in wfscent")
##             self.outputData=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,2),Numeric.Float32)       # Centroid arrays
##         else: #fpga can access whole array.
##             self.doPartialFPGA=0
##             self.fpgaarr=fpga.mallocHostMem(self.fpid,(arrsize,),Numeric.Int8)
##             if self.atmosPhaseType=="phaseonly":
##                 tmp=cmod.utils.arrayFromArray(self.fpgaarr,(self.wfs_nsubx*self.wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n),Numeric.Float32)
##                 self.reorderedPhs=Numeric.reshape(tmp,(self.wfs_nsubx,self.wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n))
##             elif self.atmosPhaseType=="phaseamp":#shouldn't get here cos error
##                 tmp=cmod.utils.arrayFromArray(self.fpgaarr,(self.wfs_nsubx*self.wfs_nsubx*self.nIntegrations,self.wfs_n,self.wfs_n,2),Numeric.Float32)
##                 self.reorderedPhs=Numeric.reshape(tmp,(self.wfs_nsubx,self.wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n,2))
##             else:
##                 raise Exception("wfscent: atmosPhaseType")
##             self.outputData=cmod.utils.arrayFromArray(self.fpgaarr[self.wfs_nsubx*self.wfs_nsubx*self.nIntegrations*self.wfs_n*self.wfs_n*4:],(self.wfs_nsubx,self.wfs_nsubx,2),Numeric.Float32)
##         self.centx=self.outputData[:,:,0]
##         self.centy=self.outputData[:,:,1]
##         #set up memory...
##         self.loadPoisson()
##         self.loadSeed()
##         self.loadPupil()

##         #now initialise the registers...
##         self.loadBinCtrl()
##         self.loadCenMask()
##         self.loadFPGAReg()
##         self.fpgaInitialised=1

##     def loadPupil(self):
##         """Load the pupil map - note, using this overwrites the reorderedPhs array..."""
##         #Now, place the pupil into a bit format that can be read by FPGA.
##         npxls=(self.wfs_nsubx*self.wfs_n)**2
##         usexsubap=self.wfs_nsubx
##         useysubap=self.wfs_nsubx
##         self.usePupil=1
##         self.symtype=0
##         if npxls>4*1024*1024*8*4:
##             print "Warning: can't fit pupil function into FPGA.  Assuming all pixels needed."
##             self.usePupil=0
##         elif npxls>4*1024*1024*8*2:
##             print "Using 2 fold symmetry for pupil function"
##             self.symtype=2
##         elif npxls>4*1024*1024*8:
##             print "Using 1 fold symmetry for pupil function"
##             self.symtype=1
##         self.puparr=Numeric.zeros((4*1024*1024,),Numeric.UInt8)
##         if self.usePupil==0:#use every pixel.
##             self.puparr[:]=0xff
##         else:
##             if self.symtype==0:
##                 pupfn=self.pupsub
##             elif self.symtype==1:
##                 useysubap=(useysubap+1)/2
##                 pupfn=self.pupsub[:useysubap]
##             else:
##                 usexsubap=(usexsubap+1)/2
##                 useysubap=(useysubap+1)/2
##                 pupfn=self.pupsub[:useysubap,:usexsubap]
##             for i in range(useysubap):
##                 for j in range(usexsubap):
##                     for k in range(self.wfs_n):
##                         for l in range(self.wfs_n):
##                             indx=l+self.wfs_n*k+self.wfs_n**2*(j+usexsubap*i)
##                             if pupfn[i,j,k,l]==1 and self.subflag[i,j]==1:
##                                 self.puparr[indx/8]=self.puparr[indx/8] | (1<<(indx%8))
##         #finished getting it in bit form... so...
##         #now check that the pupil has been created in a legal symmetrical form - if not, warn user.
##         if self.symtype==1:
##             tmppupsub=self.pupsub.copy()
##             tmppupsub[(self.wfs_nsubx+1)/2:]=self.pupsub[self.wfs_nsubx/2-1::-1]
##             if Numeric.sum((1-(tmppupsub==self.pupsub)).flat)>0:
##                 print "WARNING: pupil mask is not subap-symmetric about the x axis - the FPGA will be using a slightly different pupil mask"
##         elif self.symtype==2:
##             tmppupsub=self.pupsub.copy()
##             tmppupsub[(self.wfs_nsubx+1)/2:]=self.pupsub[self.wfs_nsubx/2-1::-1]
##             tmppupsub[:,(self.wfs_nsubx+1)/2:]=tmppupsub[:,self.wfs_nsubx/2-1::-1]
##             if Numeric.sum((1-(tmppupsub==self.pupsub)).flat)>0:
##                 print "WARNING: pupil mask is not subap-symmetric about the x axis - the FPGA will be using a slightly different pupil mask"

##         tmparr=cmod.utils.arrayFromArray(self.fpgaarr,(4*1024*1024,),Numeric.UInt8)
##         tmparr[:]=self.puparr[:]
##         #now load the data to the fpga.
##         print "Loading pupil function into FPGA..."
##         fpid=self.fpid
##         fpga.writeReg(fpid,0x2,6)#stop pipe
##         fpga.writeAddr(fpid,tmparr,1)#input address
##         fpga.writeReg(fpid,4*1024*1024/8,2)#size (quad words)
##         fpga.writeReg(fpid,0,4)#size to write (0 bytes).
##         fpga.writeReg(fpid,64,6)#reset input/output fifos.
##         fpga.writeReg(fpid,8,512)#set pipe to write pupil fn.
##         addr=fpga.readReg(fpid,525)
##         print "Loading Pupil: Current QDR address is: %s, setting to zero"%str(addr)
##         fpga.writeReg(fpid,0,525)
##         fpga.writeReg(fpid,1,6)#set it going.
##         while 1:#should only print a few messages here at most before done.
##             time.sleep(0.008)
##             addr=fpga.readReg(fpid,525)
##             print "Loading Pupil: Writing to address %s"%str(addr)
##             if addr==0 or addr==2**19:
##                 break
##         print "FPGA QDR memory filled with pupil function"
##         time.sleep(0.01)
##         fpga.writeReg(fpid,0,512)#unset the QDR write pipe.

##     def loadSeed(self):
##         """Load the random number generator seed (for readout noise).  This is taken from the Cray mta_test.c example, converted to python.  Note, this overwrites any data in reorderedPhs..."""
##         defaultSeed=4357L
##         int32_mask=0xffffffff
##         multiplier=1812433253L #Don Knuth, Vol 2
##         seedarr=cmod.utils.arrayFromArray(self.fpgaarr,(624,),Numeric.Int32)#really should be UInt32, but haven't recompiled fpgamodule.c to cope yet.  Seed different, but so what!.
##         s=defaultSeed
##         seedarr[0]=s&int32_mask
##         for i in range(1,624):
##             tmp=(multiplier*(seedarr[i-1]^(seedarr[i-1]>>30))+i)
##             seedarr[i]=tmp&int32_mask
##         #now load the data to the fpga.
##         print "Random seed memory array created, loading into FPGA..."
##         fpid=self.fpid
##         fpga.writeReg(fpid,0x2,6)#stop pipe
##         fpga.writeAddr(fpid,seedarr,1)#input address
##         fpga.writeReg(fpid,624/2,2)#size (quad words)
##         fpga.writeReg(fpid,0,4)#size to write (0 bytes).
##         fpga.writeReg(fpid,64,6)#reset input/output fifos.
##         fpga.writeReg(fpid,4,512)#set pipe to write QDR.
##         fpga.writeReg(fpid,1,6)#set it going.
##         time.sleep(0.01)
##         print "Written seed to FPGA"
##         fpga.writeReg(fpid,0,512)#unset the QDR write pipe.
        
##     def loadPoisson(self):
##         """Load the QDR memory with appropriate random variables.
##         For <2ppp, have 256 bins, each with 512 byte entries. (64 qwords)
##         For <8ppp, have 384 bins, each with 1024 byte entries. (128 qwords)
##         For <32ppp, have 768 bins, each with 2048 byte entries. (256 qwords)
##         For >=32ppp, have gaussian. (262144 qwords, 2MB, 1048576 short entries)
##           --addressing scheme: 20 bits (only 19 valid for QDR).
##           --01234567890123456789
##           --000000ppppppppcccccc     <2 (14 bits)
##           --0000aapppppppccccccc    2-8 (16 bits)
##           --00aappppppppcccccccc   8-32 (18 bits)
##           --01iiiiiiiiiiiiiiiiii   Gaussian.
##           --(aa here means that at least one of these numbers will be 1.
##           --p represents bits coming from the light level, c represents bits
##           --coming from the bin, ie the dpbm output)
##         Note, this overwrites any data in the reorderedPhs array...
##         """
##         qdrmem=cmod.utils.arrayFromArray(self.fpgaarr,(4*1024*1024,),Numeric.UInt8)#everything
##         qdrmemp=cmod.utils.arrayFromArray(self.fpgaarr,(2*1024*1024,),Numeric.UInt8)#poisson part
##         qdrmems=cmod.utils.arrayFromArray(self.fpgaarr[2*1024*1024:],(1024*1024,),Numeric.Int16)#gaussian part - signed 8.8 format.
##         for i in range(256):#less than 2ppp, fill the array.
##             ppp=i*2**-7#the mean light level of this bin.
##             #RandomArray.poisson(ppp,(512,)).astype(Numeric.UInt8)
##             if self.control["cal_source"]==0:
##                 qdrmemp[i*512:i*512+512]=self.mypoissondist(ppp,512)
##             else:
##                 qdrmemp[i*512:i*512+512]=int(ppp)
##         for i in range(384):
##             #2-8ppp.
##             ppp=2+i*2**-6
##             if self.control["cal_source"]==0:
##                 qdrmemp[256*512+i*1024:256*512+i*1024+1024]=self.mypoissondist(ppp,1024)
##             else:
##                 qdrmemp[256*512+i*1024:256*512+i*1024+1024]=int(ppp)
##         for i in range(768):
##             ppp=8+i*2**-5
##             if self.control["cal_source"]==0:
##                 qdrmemp[256*512+384*1024+i*2048:256*512+384*1024+i*2048+2048]=self.mypoissondist(ppp,2048)
##             else:
##                 qdrmemp[256*512+384*1024+i*2048:256*512+384*1024+i*2048+2048]=int(ppp)
##         if self.control["cal_source"]==0:
##             qdrmems[:]=(RandomArray.standard_normal(1024*1024)*2**8).astype(Numeric.Int16)
##         else:
##             qdrmems[:]=0
##         self.savedQDRMem=qdrmem.copy()

##         #now we have the memory done, should load it into the FPGA.
##         print "QDR memory array created, loading into FPGA..."
##         fpid=self.fpid
##         fpga.writeReg(fpid,0x2,6)#stop pipe
##         fpga.writeAddr(fpid,qdrmem,1)#input address
##         fpga.writeReg(fpid,4*1024*1024/8,2)#size (quad words)
##         fpga.writeReg(fpid,0,4)#size to write (0 bytes).
##         fpga.writeReg(fpid,64,6)#reset input/output fifos.
##         fpga.writeReg(fpid,1,512)#set pipe to write QDR.
##         addr=fpga.readReg(fpid,525)
##         print "Loading poisson RV to FPGA: Current QDR address is %s, setting to zero"%str(addr)
##         fpga.writeReg(fpid,0,525)
##         fpga.writeReg(fpid,1,6)#set it going.
##         while 1:#should only print a few messages here at most before done.
##             time.sleep(0.008)
##             addr=fpga.readReg(fpid,525)
##             print "Loading Poisson RV: Writing to address %s"%str(addr)
##             if addr==0 or addr==2**19:
##                 break
##         print "QDR memory filled"
##         fpga.writeReg(fpid,0,512)#unset the QDR write pipe.

##     def loadBinCtrl(self):
##         """Load control vectors for binning into the FPGA"""
##         #now set up the binning control...
##         tmp=Numeric.zeros((32,),Numeric.Int8)
##         cnt=0
##         #first do x...
##         for i in range(32):
##             cnt+=self.wfs_nimg
##             if cnt>=self.wfs_nfft:
##                 excess=cnt-self.wfs_nfft
##                 use=self.wfs_nimg-excess
##                 tmp[i]=(use&0x3f)|0x40
##                 cnt=excess
##             else:
##                 tmp[i]=(self.wfs_nimg&0x3f)
##         larr=cmod.utils.arrayFromArray(tmp,(4,),Numeric.Int64)
##         for i in range(4):
##             fpga.writeReg(self.fpid,larr[i],536+i)
##         cnt=0
##         #now for the y...
##         for i in range(32):
##             cnt+=self.wfs_nimg
##             if cnt>=self.wfs_nfft:
##                 excess=cnt-self.wfs_nfft
##                 use=self.wfs_nimg-excess
##                 tmp[i]=(use&0x3f)|0x40
##                 cnt=excess
##             else:
##                 tmp[i]=(self.wfs_nimg&0x3f)
##         for i in range(4):
##             fpga.writeReg(self.fpid,larr[i],540+i)
##         #have now finished setting up the binning control.
##     def loadCenMask(self):
##         """Load the centroid mask into the FPGA registers..."""
##         iarr=Numeric.zeros((32,),Numeric.UInt32)
##         larr=cmod.utils.arrayFromArray(iarr,(16,),Numeric.Int64)
##         #Largest cent mask to be loaded is 32x32 pixels.  So, have to map
##         #out cent mask onto an array this size.
##         pos=0
##         cenmask=self.cenmask.astype(Numeric.Int64)
##         for y in range(self.wfs_nimg):
##             for x in range(self.wfs_nimg):
##                 iarr[pos/32]=iarr[pos/32] | (cenmask[y,x]<<(pos%32))
##                 pos+=1
##         #print larr
##         for i in range(16):
##             fpga.writeReg(self.fpid,larr[i],544+i)
        
##     def loadFPGAReg(self):
##         """Initialise the FPGA registers"""
##         fpid=self.fpid
##         fpga.writeReg(fpid,0x2,6)#stop pipe
##         if self.doPartialFPGA:
##             fpga.writeAddr(fpid,self.fpgaInArr,1)#input address
##             fpga.writeReg(fpid,self.fittedSubaps*self.nIntegrations*self.wfs_n*self.wfs_n*4/8,2)#size (quad words)
##             fpga.writeAddr(fpid,self.fpgaOutArr,3)#output address
##             fpga.writeReg(fpid,self.fittedSubaps,4)#size to write (in bytes).
##         else:
##             fpga.writeAddr(fpid,self.reorderedPhs,1)#input address
##             fpga.writeReg(fpid,self.wfs_nsubx*self.wfs_nsubx*self.nIntegrations*self.wfs_n*self.wfs_n*4/8,2)#size (quad words)
##             fpga.writeAddr(fpid,self.outputData,3)#output address
##             fpga.writeReg(fpid,self.wfs_nsubx*self.wfs_nsubx,4)#size to write (in bytes).

##         fpga.writeReg(fpid,64,6)#reset input/output fifos.
##         fpga.writeReg(fpid,2,512)#set pipe to do WFSing.
##         fpga.writeReg(fpid,self.nIntegrations,519)#number of integrations
##         fpga.writeReg(fpid,self.wfs_nimg,520)#centroid subap x size
##         fpga.writeReg(fpid,self.wfs_nimg,521)#centroid subap y size
##         fpga.writeReg(fpid,self.wfs_nfft,522)#fft size
##         fpga.writeReg(fpid,self.wfs_n,523)#phase size
##         fpga.writeReg(fpid,self.sigFPGA,524)#user scale
##         fpga.writeReg(fpid,2**15/self.wfs_nfft,526)#inv fft size
##         #fpga.writeReg(fpid,self.subapsize,527)#scale/norm subapsize
##         fpga.writeReg(fpid,(1-self.control["cal_source"])*int(self.skybrightness*2**10/self.wfs_n**2),528)#sky brightness in 16.10 fixed point format per pupil pixel in each subap...
##         fpga.writeReg(fpid,(1-self.control["cal_source"])*int(self.wfs_floor),529)#threshold value to remove
##         fpga.writeReg(fpid,(1-self.control["cal_source"])*int((self.wfs_read_mean-self.wfs_read_sigma*127/26.11)*256),530)#mean readout noise (const level added to signal)
##         fpga.writeReg(fpid,(1-self.control["cal_source"])*int(self.wfs_read_sigma/26.11*256),531)#standard deviation of readout noise (RMS readout noise)
##         fpga.writeReg(fpid,self.wfs_nimg,532)#centroid subap x size (again)
##         fpga.writeReg(fpid,self.wfs_nimg,533)#centroid subap y size (again)
##         fpga.writeReg(fpid,(self.wfs_nsubx&0x3ff)|((self.wfs_nsubx&0x3ff)<<10),534)#write the number of subaps
##         fpga.writeReg(fpid,self.symtype,535)#write the symmetry type

##     def initDummyFPGA(self):
##         """Set up memory with same name as FPGA versions"""
##         if self.atmosPhaseType=="phaseonly":
##             self.reorderedPhs=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n),self.fpDataType)
##         elif self.atmosPhaseType=="phaseamp":
##             self.reorderedPhs=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n,2),self.fpDataType)
##         else:
##             raise Exception("AtmosPHaseType wfscent")
##         self.outputData=Numeric.zeros((self.wfs_nsubx,self.wfs_nsubx,2),Numeric.Float32)
##         # Centroid arrays
##         self.centx=self.outputData[:,:,0]
##         self.centy=self.outputData[:,:,1]
    
    def generateNext(self,msg=None):
        """WFS main iteration loop.
        Not expecting any msgs.
        """
        t1=time.time()
        if self.debug!=None:
            print "wfs: Doing generateNext (debug=%s)"%str(self.debug)
        current=self.control["parentType"]
        if self.generate==1:
            if self.newDataWaiting:#this is always 1
                if self.parent[current].dataValid==1:
                    #self.inputData=self.parent.outputData
                    self.inputInvalid=0
                else:
                    print "wfscent: waiting for data from dm, but not valid"
                    self.inputInvalid=1
                    self.dataValid=0
            if self.inputInvalid==0: # there was an input, so we can integrate...
                self.dataValid=0
                wfs=self.wfscentObj # has been copied from thisObjList before generateNext is called...
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
                    print "wfs time:",time.time()-t1
            #self.config.setVal("simulationTime",self.config.getVal("simulationTime")+self.tstep)
        else:
            self.dataValid=0
        if self.debug!=None:
            print "wfs: Done generateNext (debug=%s)"%str(self.debug)
        self.generateNextTime=time.time()-t1








##     def runFPGA(self):
##         """Tell the FPGA where the data is..."""
##         fpid=self.fpid
##         if self.control["cal_source"]!=self.lastCalSource:
##             #user has requested a change since last iteration - have to reload poisson/readnoise.
##             self.loadPoisson()
##             fpga.writeReg(fpid,(1-self.control["cal_source"])*int((self.wfs_read_mean-self.wfs_read_sigma*127/26.11)*256),530)#mean readout noise (const level added to signal)
##             fpga.writeReg(fpid,(1-self.control["cal_source"])*int(self.wfs_read_sigma/26.11*256),531)#standard deviation of readout noise (RMS readout noise)
##             fpga.writeReg(fpid,(1-self.control["cal_source"])*int(self.skybrightness*2**10/self.wfs_n**2),528)
##             fpga.writeReg(fpid,(1-self.control["cal_source"])*int(self.wfs_floor),529)
##             #now reset the read/write DMA addresses...
##             fpga.writeReg(fpid,0x2,6)#stop pipe
##             if self.doPartialFPGA:
##                 fpga.writeAddr(fpid,self.fpgaInArr,1)#input address
##                 fpga.writeReg(fpid,self.fittedSubaps*self.nIntegrations*self.wfs_n*self.wfs_n*4/8,2)#size (quad words)
##                 fpga.writeAddr(fpid,self.fpgaOutArr,3)#output address
##                 fpga.writeReg(fpid,self.fittedSubaps,4)#size to write (in bytes).
##             else:
##                 fpga.writeAddr(fpid,self.reorderedPhs,1)#input address
##                 fpga.writeReg(fpid,self.wfs_nsubx*self.wfs_nsubx*self.nIntegrations*self.wfs_n*self.wfs_n*4/8,2)#size (quad words)
##                 fpga.writeAddr(fpid,self.outputData,3)#output address
##                 fpga.writeReg(fpid,self.wfs_nsubx*self.wfs_nsubx,4)#size to write (in bytes).
##             fpga.writeReg(fpid,64,6)#reset input/output fifos (prob not needed).
##             fpga.writeReg(fpid,2,512)#set pipe to do WFSing.
                
##         self.lastCalSource=self.control["cal_source"]

##         #now, we set the FPGA going (after copying data if necessary).
##         if self.doPartialFPGA:
##             #array is too large to DMA all at once to FPGA, so do in parts.
##             t0=time.time()
##             #calculate number of times a full array is needed...
##             partialFull=int(self.wfs_nsubx*self.wfs_nsubx/self.fittedSubaps)
##             #and the number of subaps left to do.
##             leftOver=self.wfs_nsubx*self.wfs_nsubx-partialFull*self.fittedSubaps
##             waitLeftOver=leftOver*self.nIntegrations*self.wfs_n*self.wfs_n*5e-9
##             if self.atmosPhaseType=="phaseonly":
##                 reordered=cmod.utils.arrayFromArray(self.reorderedPhs,(self.wfs_nsubx*wfs_nsubx,self.nIntegrations,self.wfs_n,self.wfs_n),Numeric.Float32)
##             else:
##                 raise Exception("not phaseonly")
##             for i in range(partialFull):
##                 #copy memory into FPGA buffer
##                 self.fpgaInArr[:]=reordered[i*fittedSubaps:(i+1)*fittedSubaps]
##                 fpga.writeReg(fpid,0x2,6)#reinitialise
##                 fpga.writeReg(fpid,1,6)#set it going.
##                 time.sleep(self.waitFPGATime/partialFull)#wait for it to complete (or almost)
##                 while fpga.readReg(fpid,5)!=7:#wait for reading to complete by checking register.
##                     pass
##             if leftOver>0:#copy the last bit...
##                 self.fpgaInArr[:leftOver]=reordered[partialFull*fittedSubaps:partialFull*fittedSubaps+leftOver]
##                 fpga.writeReg(fpid,leftOver*self.nIntegrations*self.wfs_n*self.wfs_n*4/8,2)#read siz
##                 fpga.writeReg(fpid,leftOver,4)#size to write
##                 fpga.writeReg(fpid,0x2,6)#reinitialise
##                 fpga.writeReg(fpid,1,6)#set it going
##                 time.sleep(waitLeftOver)
##                 while fpga.readReg(fpid,5)!=7:#wait til finished
##                     pass
##                 #now reset the registers for next time...
##                 fpga.writeReg(fpid,self.fittedSubaps*self.nIntegrations*self.wfs_n*self.wfs_n*4/8,2)#read siz
##                 fpga.writeReg(fpid,self.fittedSubaps,4)#size to write

##             t1=time.time()
##             if self.timing:
##                 print "WFSCent time taken: %s"%str(t1-t0)
                      
##         else:
##             #all subaps at once...
##             #fpga.writeReg(fpid,0x2,6)#stop the pipe (this may not be necessary).
##             t0=time.time()
##             #self.outputData[-1,-1,0]=1.
##             #tmp=float(self.outputData[-1,-1,0])
##             fpga.writeReg(fpid,0x2,6)#reinitialise
##             #ttmp=[hex(fpga.readReg(self.fpid,5))]
##             fpga.writeReg(fpid,1,6)#set it going.
##             #ttmp.append(hex(fpga.readReg(self.fpid,5)))
##             if self.waitFPGA:
##                 time.sleep(self.waitFPGATime)#wait for it to complete (or almost).
##                 while fpga.readReg(fpid,5)!=7:#wait for reading to complete by checking register...
##                     pass
##                     #ttmp.append(hex(fpga.readReg(self.fpid,5)))
## #            while tmp==self.outputData[-1,-1,0]:
## #                pass
##             t1=time.time()
##             if self.timing:
##                 print "WFSCent Time taken: %s"%str(t1-t0)



            
##     def reorder(self,phs,pos):
##         """Do a reodering of the phase buffer, so that it is in the form ready
##         for the fpga.  Phases corresponding to subaps should be placed next to
##         those for the same subap at a later time (greater pos).
##         Also use with createSHImg() but not integrate().  If the FPGA has been initialised, this will place it in FPGA readable memory..."""
##         #nsubx=    self.wfs_nsubx
##         n=    self.wfs_n
##         typecode=self.reorderedPhs.typecode()
##         if self.atmosPhaseType=="phaseonly":
##             for i in range(self.wfs_nsubx):
##                 for j in range(self.wfs_nsubx):
##                     # Subap phase array
##                     self.reorderedPhs[i,j,pos,:,:]=phs[i*n:(i+1)*n,j*n:(j+1)*n].astype(typecode)
##         elif self.atmosPhaseType=="phaseamp":
##             for i in range(self.wfs_nsubx):
##                 for j in range(self.wfs_nsubx):
##                     # Subap phase array
##                     self.reorderedPhs[i,j,pos,:,:,0]=phs[0,i*n:(i+1)*n,j*n:(j+1)*n].astype(typecode)#phase
##                     self.reorderedPhs[i,j,pos,:,:,1]=phs[1,i*n:(i+1)*n,j*n:(j+1)*n].astype(typecode)#amplitude
##         else:
##             raise Exception("atmosphasetype")

##     def createSHImgs(self):
##         """Do the FFTs etc and add integrations to create the powerspectra (high light level images).
##         All this will eventually (we hope) be done in the FPGAs"""
##         tmp=0.5*float(self.wfs_n)/float(self.wfs_nfft)*2.*math.pi
##         self.subimg*=0.0                                                  # Zero the CCD
##         #nsubx=self.wfs_nsubx
##         for i in range(self.wfs_nsubx):
##             for j in range(self.wfs_nsubx):
##                 if self.subflag[i][j]==1:
##                     for k in range(self.nIntegrations):
##                         if self.atmosPhaseType=="phaseonly":
##                             phssub=(self.reorderedPhs[i,j,k]-tmp*self.xtiltfn-tmp*self.ytiltfn).astype(self.fpDataType)
##                         elif self.atmosPhaseType=="phaseamp":
##                             phssub=(self.reorderedPhs[i,j,k,:,:,0]-tmp*self.xtiltfn-tmp*self.ytiltfn).astype(self.fpDataType)
##                         #now do the FFT: plan, phsin, imgout, pypupil
##                         if self.fpDataType==Numeric.Float64:
##                             if phssub.typecode()!=self.fpDataType or self.tmpImg.typecode()!=self.fpDataType or self.pupsub.typecode()!=self.fpDataType:
##                                 print "ERROR with typecodes in wfs.py"
##                                 raise Exception("Error with typecodes in wfs.py mkimg")
##                             if self.atmosPhaseType=="phaseonly":
##                                 cmod.mkimg.mkimg(self.fftwPlan,phssub,self.tmpImg,self.pupsub[i][j])
##                             elif self.atmosPhaseType=="phaseamp":
##                                 cmod.mkimg.mkimg(self.fftwPlan,phssub,self.tmpImg,(self.pupsub[i,j]*self.reorderedPhs[i,j,k,:,:,1]).astype(self.fpDataType))
##                             elif self.atmosPhaseType=="realimag":
##                                 raise Exception("realimag")
##                         else:
##                             if phssub.typecode()!=self.fpDataType or self.tmpImg.typecode()!=self.fpDataType or self.pupsub.typecode()!=self.fpDataType:
##                                 print "ERROR with typecodes in wfs.py",phssub.typecode(),self.tmpImg.typecode(),self.pupsub.typecode()
##                                 raise Exception("Error with typecodes in wfs.py mkimg")
##                             if self.atmosPhaseType=="phaseonly":
##                                 cmod.mkimgfloat.mkimg(self.fftwPlan,phssub,self.tmpImg,self.pupsub[i][j])
##                             elif self.atmosPhaseType=="phaseamp":
##                                 cmod.mkimgfloat.mkimg(self.fftwPlan,phssub,self.tmpImg,(self.pupsub[i,j]*self.reorderedPhs[i,j,k,:,:,1]).astype(self.fpDataType))
##                             elif self.atmosPhaseType=="realimag":
##                                 raise Exception("realimag")
                            
##                         self.subimg[i][j]+=self.tmpImg                    # Long exposure image




##     def integrate(self,phs):
        
##         """Shack Hartmann image integration.
##         Use either this OR reorder and createSHImg"""
##         n=    self.wfs_n
##         #nfft=    self.wfs_nfft
##         #nsubx=    self.wfs_nsubx
##         tmp=0.5*float(n)/float(self.wfs_nfft)*2.*math.pi
##         #tmp=0.
##         #print "INTEGRATING",self.texp
##         phssub=Numeric.zeros((n,n),self.fpDataType)#Numeric.Float64
##         #agb: This is probably where we can start putting into FPGA
##         #Integration may be a problem - can we store partially
##         #integrated subaps in the SRAM or something?  How big would
##         #they need to be (binning could probably be done before
##         #integrating if needed).
        
##         for i in range(self.wfs_nsubx):                                                                              # Loop over the subapertures
##             for j in range(self.wfs_nsubx):
##                 if(self.subflag[i][j]==1):        
##                     phssub=(phs[i*n:(i+1)*n,j*n:(j+1)*n].astype(self.fpDataType)).astype(self.fpDataType)   # Subaperture phase array
##                     #if(sum(sum(phssub))!=0.):# Speed up for poking - don't make image 
##                     phssub=phssub-tmp*self.xtiltfn-tmp*self.ytiltfn                                         # Centers image on 4 pixels if phase=0
                    
##                     #now do the FFT
##                     if self.fpDataType==Numeric.Float64:
##                         if phssub.typecode()!=self.fpDataType or self.tmpImg.typecode()!=self.fpDataType or self.pupsub.typecode()!=self.fpDataType:
##                             print "ERROR with typecodes in wfs.py"
##                             raise Exception("Error with typecodes in wfs.py mkimg")
##                         cmod.mkimg.mkimg(self.fftwPlan,phssub,self.tmpImg,self.pupsub[i][j])                # phsin, imgout, pypupil
##                     else:
##                         if phssub.typecode()!=self.fpDataType or self.tmpImg.typecode()!=self.fpDataType or self.pupsub.typecode()!=self.fpDataType:
##                             print "ERROR with typecodes in wfs.py"
##                             raise Exception("Error with typecodes in wfs.py mkimg")
##                         cmod.mkimgfloat.mkimg(self.fftwPlan,phssub,self.tmpImg,self.pupsub[i][j])           # phsin, imgout, pypupil
                        
##                     self.subimg[i][j]+=self.tmpImg


                    
                    
##     def tidyImage(self):
##         """Flip image (to correct for FFT), then bin the image up if
##         oversampled, then add read and shot noise."""
        
##         self.shimg*=0.                                                # Reset...=zeros((nsubx*nimg,nsubx*nimg),Float)
##         readnoise=  self.wfs_read_sigma
##         mean=self.wfs_read_mean
## ##         m=0
## ##         mn=1000000
##         #nsubx=self.wfs_nsubx
##         for i in range(self.wfs_nsubx):                                        # Loop over subaps
##             for j in range(self.wfs_nsubx):
##                 if(self.subflag[i][j]==1):
##                     bimg=self.bimg[i][j]
##                     img=fliparray(self.subimg[i][j])                          # Get subap image, flip it
##                     #img=self.conv(img)                                       # Convolve with a PSF  (include LGS psf here)
##                     cmod.binimg.binimg(img,bimg)                              # Bin it up
##                     totsig = Numeric.sum(Numeric.sum(bimg))
##                     nphs=float(self.subarea[i,j])/self.wfs_n**2#fraction of active phase pixels
##                     if(totsig>0.):
##                         bimg*=self.sig*nphs/totsig  #rescale the image.
##                         #Note, used to be an error here - should also scale by the number of pixels that receive photons.  Fixed...
##                         if self.atmosPhaseType!="phaseonly":
##                             print "Scaling of SH image? Is it needed when have phase an amplitude for atmosphere"
##                     #now add sky brightness scaled by the number of phase pixels that are in the aperture:
##                     if not self.control["cal_source"]:                      
##                         bimg+=self.skybrightness*nphs
##                         if totsig>0. or self.skybrightness*nphs>0.:
##                             # Inject shot noise
##                             cmod.imgnoise.shot(bimg,bimg)
##                         # Generate random read noise :
##                         bimg+=RandomArray.normal(mean,readnoise,bimg.shape)

##                     self.shimg[i*self.wfs_nimg:(i+1)*self.wfs_nimg,j*self.wfs_nimg:(j+1)*self.wfs_nimg]=bimg    # Tessalate up for WFS display

##     def calc_cents(self):
##         """Centroid calculation:
##         Subtracts noise background
##         Computes centroids
##         No return value."""
##         nfft=   self.wfs_nfft
##         nimg=self.wfs_nimg
##         nsubx=  self.wfs_nsubx
##         floor=  self.wfs_floor
##         #read=  self.wfs_read
##         indx=   self.tilt_indx
##         #bimg=self.bimg
##         self.outputData[:]=0.                                      # Reset...
##         #self.centx=zeros((nsubx,nsubx),Float)
##         #self.centy=zeros((nsubx,nsubx),Float)
##         #self.shimg*=0.#reset...=zeros((nsubx*nimg,nsubx*nimg),Float)
        
##         for i in range(nsubx):                                  # Loop over subaps
##             for j in range(nsubx):
##                 if(self.subflag[i][j]==1):
##                     #threadhold the SH images.
##                     if not self.control["cal_source"]:
##                         cimg=Numeric.where(self.bimg[i][j]<self.wfs_floor,0,self.bimg[i][j]-self.wfs_floor)*self.cenmask
##                     else:
##                         cimg=self.bimg[i][j]*self.cenmask
##                     #bimg=self.bimg[i][j]                        # This will be only integers...
##                     #img=fliparray(self.subimg[i][j])# Get subap image
##                     #img=self.conv(img)# Convolve with a PSF
##                     #binimg.binimg(img,bimg)# bin it up
##                     #totsig = Numeric.sum(bimg.flat) 
##                     #if(totsig>0):
##                     #if(not(self.control['cal_source'])):#agb why?
##                     #    img0=((bimg/totsig)*self.sig)+read*read
##                     #    imgnoise.shot(img0,bimg)# Inject shot noise
##                     #bimg=Numeric.clip(bimg,floor,1.e6)      # Floor==value in electrons
##                     #cimg=(bimg-floor)*self.cenmask          # Apply threshld,centroding mask
##                     totsig = Numeric.sum(Numeric.sum(cimg))
##                     if(totsig==0.):
##                         totsig=1.#division by zero...
##                     # Centroid calculation
##                     self.centx[i,j]=Numeric.sum(Numeric.sum(cimg)*indx)/totsig  
##                     self.centy[i,j]=Numeric.sum(Numeric.sum(Numeric.transpose(cimg))*indx)/totsig
                        

##     def tilt(self,x,y):
##         return y

    
##     def wfs_sig(self):
##         """wfs signal per exposure"""
##         bandwidth=self.wfs_bandwidth                                  # Optical bandwidth (Angstrom)
##         thruput=self.wfs_thruput                                      # Optical thruput incl. DQE
##         rate_factor=self.wfs_phot_rate_factor                         # Photons/cm^2/s/Angstrom
##         pupil_area=(self.telDiam/float(self.wfs_nsubx)*100.)**2.      # Pupil area/cm^2
##         sig=rate_factor*thruput*bandwidth*self.wfs_int*pupil_area/(2.5**self.wfs_mag) # Detected image photons
##         #print 'WFS  photons/subap/integration: ',sig
##         return sig

    
##     def conv(self,img):
##         """Convolve subap image with a PSF
##         e.g. LGS spot shape"""
##         nfft=self.wfs_nfft
##         temp=Numeric.zeros((nfft,nfft),Numeric.Float64)
##         temp[nfft/2-1,nfft/2-1]=1.
##         temp[nfft/2,nfft/2]=1.
##         temp[nfft/2-1,nfft/2]=1.
##         temp[nfft/2,nfft/2-1]=1.
##         wfsimg=fliparray(temp)
##         temp1=FFT.real_fft2d(img)
##         temp2=FFT.real_fft2d(wfsimg)
##         convimg=FFT.inverse_real_fft2d(temp1*temp2)
##         return convimg

##     def mypoissondist(self,mean,nsamp=1024):
##         """Generates a poisson distribution by considering the probability distribution and then randomising the order.  This ensures that even after the sample has been used many times, the overall probability is still correct.  See py/poisson/poisson.py for more details."""
##         if mean>32:
##             raise Exception("Mean value too large, use normal distribution")
##         y=[]
##         x=range(64)
##         out=[]
##         for i in range(64):
##             tmp=int(round(self.calcpprob(i,mean)*nsamp))
##             y.append(tmp)
##             for j in range(tmp):
##                 out.append(i)
##         #print len(out)
##         while len(out)<nsamp:
##             out.append(int(round(mean)))
##         out=out[:nsamp]
##         #now randomise the order...
##         d={}
##         for i in out:
##             k=random.random()
##             while d.has_key(k):
##                 k=random.random()
##             d[k]=i
##         l=d.keys()
##         l.sort()
##         out=[]
##         for i in l:
##             out.append(d[i])
##         return Numeric.array(out,Numeric.UnsignedInt8)
##     def fact(self,x):
##         f=1.
##         i=1.
##         while i<=x:
##             f*=i
##             i+=1
##         return f
##     def calcpprob(self,x,mean):
##         if mean<=0 or x<0:
##             return 0
##         return math.exp(-mean)*(float(mean)**x)/self.fact(x)

    def drawCents(self,fromcent=0,objNumber=None,mask=None):
        """Draw centroids in a format that can be easily displayed.
        If fromcent is 1, will construct an image from centroid locations,
        otherwise, will use the noisy image (non-fpga only).
        objNumber specifies which resource sharing object to use.  None means
        the one currently in use.
        If mask is set, parts of the images that aren't used in centroid computation will be masked out (ncen).
        If mask is 1, they are masked at max value, and if -1, masked at zero.
        """
        if mask==None:
            mask=self.imgmask
        if objNumber==None:
            try:
                wfsobj=self.wfscentObj
            except:
                print "ERROR getting wfscentObj - assuming  thisObjList[0].wfscentObj"
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
        if mask==None:
            mask=self.imgmask
        if objno==None:
            wfsobj=self.wfscentObj
        else:
            wfsobj=self.thisObjList[objno].wfscentObj
        if img=="corrimg":
            img=wfsobj.corrimg
        elif img=="corrPattern":
            img=wfsobj.corrPatternUser
        for i in xrange(wfsobj.nsubx):
            for j in xrange(wfsobj.nsubx):
                 # Tessalate up for WFS display:
                self.shimg[i*wfsobj.nimg:(i+1)*wfsobj.nimg,j*wfsobj.nimg:(j+1)*wfsobj.nimg]=img[i,j]
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
            txt+="""<plot title="WFS SH img%s" cmd="data=%s.drawCents(0)" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            if self.imageOnly==0:
                txt+="""<plot title="XCentroids%s" cmd="data=%s.wfscentObj.outputData[:,:,0]" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
                txt+="""<plot title="YCentroids%s" cmd="data=%s.wfscentObj.outputData[:,:,1]" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
                txt+="""<plot title="1D centroids%s" cmd="data=%s.wfscentObj.outputData.ravel()" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            elif self.imageOnly==1:
                txt+="""<plot title="Centroids%s" cmd="data=%s.wfscentObj.outputData.ravel()" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            else:
                txt+="""<plot title="Centroids%s" cmd="data=%s.wfscentObj.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
            txt+="""<plot title="Centroids 2D%s" cmd="data=%s.drawCents(1)" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            for i in range(len(self.thisObjList)):
                this=self.thisObjList[i]
                wfs=this.wfscentObj
                if type(this.laserGuideStar)!=type(None):
                    if type(this.laserGuideStar)==numpy.ndarray:
                        txt+="""<plot title="LGS elongation%s" cmd="data=%s.thisObjList[%d].laserGuideStar" ret="data" type="pylab" when="cmd" palette="gray"/>\n"""%(id,objname,i)
                    else:
                        txt+="""<plot title="LGS elongation%s" cmd="data=%s.thisObjList[%d].laserGuideStar.subapImage" ret="data" type="pylab" when="cmd" palette="gray"/>\n"""%(id,objname,i)
            for i in range(len(self.thisObjList)):
                this=self.thisObjList[i]
                wfs=this.wfscentObj
                if wfs.correlationCentroiding:
                    txt+="""<plot title="Correlation%s" cmd="data=%s.drawCorrelation('corrimg')" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
                    txt+="""<plot title="Correlation PSF%s" cmd="data=%s.drawCorrelation('corrPattern')" ret="data" type="pylab" when="cmd" palette="gray"/>"""%(id,objname)

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
        #paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="dataType",typ="eval",val="this.globals.fpDataType",comment="Array numpy data type"))
        paramList.append(base.dataType.dataType(description="timing",typ="i",val="0",comment="Timing information"))
        #paramList.append(base.dataType.dataType(description="degRad",typ="eval",val="2*numpy.pi/360.",comment="degrees to radians."))        
        #paramList.append(base.dataType.dataType(description="arcsecRad",typ="eval",val="2*numpy.pi/360/3600",comment="arcsec to radians."))        
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="ntel",typ="eval",val="this.globals.npup",comment="Pixels for telescope"))
        paramList.append(base.dataType.dataType(description="telSec",typ="f",val="8.",comment="TODO: Telescope secondary diameter (m)"))
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam,wfs_nsubx,wfs_minarea)",comment="Telescope pupil"))
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))
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
        return paramList


