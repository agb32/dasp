#$Id: science.py,v 1.77 2012/01/20 06:38:23 ali Exp $
"""Module science.py : computes AO corrected PSF and associated parameters like SR, light coupling in circular and square apertures, etc...
"""

### Science camera simulation class ######################################

### Integrate AO corrected PSF and measure various parameters
#import base.readConfig#,cmod.utils
import base.aobase as aobase
#import sys,math,thread
#import math#,FFT
import os
import numpy
import string
#import cPickle
#import base.fwdMsg
#from sys import argv

#from util.dist import dist

##from cmod.psfparams import azav

#import types
import time

#from plwf import dif_
## from gist import *
## window(1,wait=1,dpi=100)
## animate(1)

##The next packages are used to write FITS files
#import util.pyfits as pyfits
import util.FITS
import util.sci
#import os,numarray


##for plotting
##animate(1)


## for debug purpose : creation of a function to compute psf from phase with numeric FFT
def phi2psf(phi,pup):
    """ Create a psf from phi and pup arrays
        Assumes Nyquist sampling"""
    import numpy.fft
    from util.flip import fliparray2 ##fliparray2 is consistent with FFT coordinate definition
    dimx=phi.shape[0] ##size of the input arrays (assumes square dimensions)
    amp=numpy.zeros((2*dimx,2*dimx),dtype="D") ## complex amplitude in the pupil

    ##we fill the complex amplitute array
    amp.real[:dimx,:dimx]=pup*numpy.cos(phi)
    amp.imag[:dimx,:dimx]=pup*numpy.sin(phi)

    ## we compute the complex amplitude in the focal plane
    amp_foyer=numpy.fft.fft2(amp)
    psf=(numpy.absolute(amp_foyer))**2
    psf=fliparray2(psf)

    ## we normalise the flux
    psf/=numpy.sum(numpy.sum(psf))

    return psf




class science(aobase.aobase):
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Create a science object which collects info such as psf, strehl etc.
        There are two parameters, science_integrate and zero_science which
        can be changed using the gui to determine whether or not to
        integrate/zero the psf signal.  These were part of the controldict
        in the old simulation code.
        Amongst other things, parent can be atmos, or dm.

        An atmosPhaseType variable can be specified in the parameter
        file, which should be one of phaseonly, phaseamp, realimag.
        phaseonly means just the atmospheric phase is passed.
        phaseamp means that the atmospheric phase and amplitude is
        passed (phs[0]==phase, phs[1]==amp array, and realimag means
        the wavefront is passed as real and imaginary parts.  In this
        case, the array should be complex.

        """
        aobase.aobase.__init__(self,parent,config,args=args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        ##Extraction of XML parameter file
        #self.debug=debug
        #self.dataValid=0
        #self.newDataWaiting=1
        #self.outputData=None#no output data...
        #self.generate=1
        #self.source=0 ##to be changed. Souce should be an element of the constructor
#        self.atmosPhaseType=config.getVal("atmosPhaseType",default="phaseonly")
#        if self.atmosPhaseType not in ["phaseonly","phaseamp","realimag"]:
#            raise Exception("science: atmosPhaseType not known %s"%self.atmosPhaseType)
#        self.sci_lam=config.getVal("sourceLam")##Imaging wavelength
        #self.plots=config.getVal("plots")
        #self.objs=config.getVal("objs") ## not really well defined
        self.timing=config.getVal("timing",default=0) ##for debugging purpose
        self.doneFinalInit=0
        self.nimgNotNfft=0
        self.ended=1
#        npup=self.npup=config.getVal("npup") ##number of linear pixels used to simulate  the pupil
#        nfft=self.nfft=self.sci_nfft=config.getVal("scinfft") ##linear number of pixels used to simulate the PSF
#        self.pupil=config.getVal("pupil") ##telescope entrance pupil
#        self.dmpupil=config.getVal("dmpupil") ##DM pupil
#        self.pup=self.pupil.fn*self.dmpupil# Pupil seen by Science camera
        #self.pup=self.pup.astype("d") ##we convert it into Float64 to be compatible with mkimg C module
#        self.pup=self.pup.astype(Numeric.Int8)
#        self.pupsum=Numeric.sum(Numeric.sum(self.pup))        
#        self.idxPup=Numeric.nonzero(self.pup.flat) #Flat version of the pixels belonging to the pupil#
#        self.telDiam=config.getVal("telDiam") ##Telescope diameter 
#        self.asRad=config.getVal("arcsecRad") ##Conversion rad-> arcsec
#        self.L_D= ( (self.sci_lam*1.e-9)/self.telDiam )/self.asRad ##Diffraction limited resolution
#        self.pix_scale=self.L_D*float(self.npup)/float(nfft) ##pixel scale (arcsec/pixel in the science image)

        self.control={}
        #leave control as global (shared between all resource sharing objects)
        self.control["zero_science"]=config.getVal("zero_science",default=10) ##idem - zero for 10 iterations...
        self.control["science_integrate"]=config.getVal("science_integrate",default=1)
        self.control["calcRMS"]=self.config.getVal("science_calcRMS",default=0)##idem
        self.control["viewCentral"]=1#if set, the gui will only be sent the central 20% of the psfs.
        self.control["lucky_integrate"]=config.getVal("lucky_integrate",default=0)
        self.thisiter=0
        self.sentPlotsCnt=0

        #self.initFPGA=self.config.getVal("initFPGA",default=0)#whether to initialise the FPGA... (load binary, set up arrays etc).
        useFPGA=self.config.getVal("useFPGA",default=0)#whether to use the FPGA initially (global)
        if useFPGA:
            self.FPGABitFile=self.config.getVal("FPGASciBitFile",default=string.join(__file__.split("/")[:-2]+["fpga","scifft.bin.ufp"],"/"))
            self.ignoreFPGAOpenFailure=self.config.getVal("ignoreFPGAOpenFailure",default=0)
        self.control["useFPGA"]=useFPGA
        self.fpDataType=self.config.getVal("fpDataType",default="f")
        self.cpDataType=self.fpDataType.upper()

        self.nthreads=self.config.getVal("nthreads",default="all")#usually an integer... or "all"
        if self.nthreads=="all":#use all available CPUs...
            self.nthreads=self.config.getVal("ncpu")#getCpus()
            print "science: Using %d threads"%self.nthreads

        for i in xrange(len(self.idstr)):
            idstr=self.idstr[i]
            parent=self.parentList[i]
            self.initialise(parent,idstr)


    def initialise(self,parent,idstr):
        this=aobase.resourceSharer(parent,self.config,idstr,self.moduleName)
        self.thisObjList.append(this)
        apt=this.config.getVal("atmosPhaseType",default="phaseonly")
        tstep=this.config.getVal("tstep")
        npup=this.config.getVal("npup") ##number of linear pixels for pupil
        nfft=this.config.getVal("scinfft",default=int(2**numpy.ceil(numpy.log2(2*npup)))) ##linear number of pixels for PSF
        nimg=this.config.getVal("scinimg",default=nfft)
        if nimg!=nfft:
            self.nimgNotNfft=1
        pupil=this.config.getVal("pupil") ##telescope entrance pupil
        realPupil=this.config.getVal("realPupil",default=None,raiseerror=0)
        sciPath=this.config.getVal("sciPath",default=idstr)
        if type(realPupil)!=type(None):
            realPupil=realPupil.fn
        try:
            dmpupil=pupil.dmpupil
        except:
            dmpupil=this.config.getVal("dmpupil",raiseerror=0,default=None) ##DM pupil
        print "todo: science - is it correct to use dmpupil here - should we not just use pupil?  dmpupil would surely be unphysical? (ask FA)"
        usedmpup=this.config.getVal("usedmpup",default=1)
        if usedmpup:
            pup=numpy.array((pupil.fn*dmpupil),numpy.int8)#.astype(numpy.int8)# Pupil seen by Science cam
        else:
            pup=pupil.fn
        #print type(pup),type(pupil.fn),type(dmpupil),pup.typecode(),pupil.fn.typecode(),dmpupil.typecode()
        atmosGeom=this.config.getVal("atmosGeom",default=None,raiseerror=0)
        sci_lam=None
        phsLam=None
        if atmosGeom!=None:
            sci_lam=atmosGeom.sourceLambda(sciPath)
            phsLam=atmosGeom.phaseLambda(sciPath)
        if sci_lam==None:
            sci_lam=this.config.getVal("sourceLam")##Imaging wavelength
            print "Warning - sci_lam not in atmosGeom - using %g"%sci_lam
        if phsLam==None:
            phsLam=sci_lam
        print "science - Using wavelength of %g nm (phase at %g nm)"%(sci_lam,phsLam)
        phaseMultiplier=phsLam/sci_lam
        telDiam=this.config.getVal("telDiam") ##Telescope diameter 
        asRad=numpy.pi/180./3600.#this.config.getVal("arcsecRad") ##Conversion rad-> arcsec
        L_D= ( (sci_lam*1.e-9)/telDiam )/asRad ##Diffraction limited resolution
        pix_scale=L_D*float(npup)/float(nfft)*float(nfft)/float(nimg) ##pixel scale (arcsec/pixel in the science image).  
        ##name of the FITS file to store the PSF
        fitsFilename=this.config.getVal("scifitsFilename",default=None,raiseerror=0)
        sciFilename=this.config.getVal("scicsvFilename",default=None,raiseerror=0)
        sciPSFSamp=this.config.getVal("sciPSFSamp",default=1)
        scinSamp=this.config.getVal("scinSamp")
        luckyFilename=this.config.getVal("luckyFilename",default=None,raiseerror=0)
        luckyImgFilename=this.config.getVal("luckyImgFilename",default=None,raiseerror=0)
        luckyImgSize=this.config.getVal("luckyImgSize",default=None,raiseerror=0)
        luckyNSampFrames=this.config.getVal("luckyNSampFrames",default=1)
        luckyHistorySize=this.config.getVal("luckyHistorySize",default=None,raiseerror=0)
        if luckyHistorySize==None:
            luckyHistorySize=int(this.config.getVal("AOExpTime")/(tstep*sciPSFSamp*luckyNSampFrames))
        luckyByteswap=this.config.getVal("luckyByteswap",default=0)
        saveFileString=this.config.getVal("scisaveFileString",raiseerror=0)#an optional string that can be used to identify a given simulation in the saved results... 
        diffPsfFilename=this.config.getVal("sciDiffPsfFilename",default=None,raiseerror=0)
        keepDiffPsf=this.config.getVal("keepDiffPsf",default=0)#overwrite mem?
        scienceListsSize=this.config.getVal("hist_list_size")
        inboxDiamList=this.config.getVal("inboxDiamList",default=[0.2])
        #check FPGA stuff...
        useFPGA=this.config.getVal("useFPGA",default=0)#whether to use the FPGA initially (global)
        if useFPGA:
            waitFPGA=this.config.getVal("waitFPGA",1)
            waitFPGATime=this.config.getVal("waitFPGATime",default=1e-6)#nfft*nfft*5e-9)
        else:
            waitFPGA=1
            waitFPGATime=1e-6
        histFilename=this.config.getVal("histFilename",default=None,raiseerror=0)
        #now create the science object that will do most of the work...
        this.sciObj=util.sci.science(npup,nfft,pup,nimg=nimg,atmosPhaseType=apt,tstep=tstep,keepDiffPsf=keepDiffPsf,pix_scale=pix_scale,fitsFilename=fitsFilename,diffPsfFilename=diffPsfFilename,scinSamp=scinSamp,sciPSFSamp=sciPSFSamp,scienceListsSize=scienceListsSize,debug=self.debug,timing=self.timing,allocateMem=0,realPup=realPupil,useFPGA=useFPGA,waitFPGA=waitFPGA,waitFPGATime=waitFPGATime,fpDataType=self.fpDataType,inboxDiamList=inboxDiamList,sciFilename=sciFilename,saveFileString=saveFileString,nthreads=self.nthreads,histFilename=histFilename,phaseMultiplier=phaseMultiplier,luckyNSampFrames=luckyNSampFrames,luckyFilename=luckyFilename,luckyImgFilename=luckyImgFilename,luckyImgSize=luckyImgSize,luckyHistorySize=luckyHistorySize,luckyByteswap=luckyByteswap)

        

        
    def finalInitialisation(self):
        """This gets called just before the main loop is entered - to finalise setting up of arrays.  We now know how large they should be, since all resource shares have been added.
        """
        #First find the max fft and pupil size...
        if self.doneFinalInit:
            return
        self.doneFinalInit=1
        nfftmax=0
        npupmax=0
        pupmult=1
        fpgarequired=0
        nimgmax=0
        fpgaObj=None
        fpgaarr=None
        for this in self.thisObjList:
            if this.sciObj.nfft>nfftmax:
                nfftmax=this.sciObj.nfft
            if this.sciObj.npup>npupmax:
                npupmax=this.sciObj.npup
            if this.sciObj.atmosPhaseType in ["phaseamp","realimag"]:
                pupmult=2
            if this.sciObj.canUseFPGA==1:
                fpgarequired=1
                fpgaObj=this.sciObj
            if this.sciObj.nimg>nimgmax:
                nimgmax=this.sciObj.nimg
        npup=npupmax
        nfft=nfftmax
        nimg=nimgmax
        if fpgarequired:
            fpid,fpgaInfo=fpgaObj.initialiseFPGA(ignoreFailure=self.ignoreFPGAOpenFailure,fpgaBitFile=self.FPGABitFile)#load fpga binary...
            fpgaarr=fpgaObj.setupFPGAArray(fpid,nfftmax)
        else:
            fpid=fpgaInfo=fpgaarr=None

        
        #now allocate memory at the max size needed.  These memory arrays will then be shared between all the sciObj objects.

        self.pupilAmplitudeMem=numpy.zeros((nfft,nfft),numpy.complex64)#self.cpDataType) ##Complex amplitude in the pupil - was complex64
        self.focusAmplitudeMem=numpy.zeros((nfft,nfft),numpy.complex64)#self.cpDataType) ##Complex amplitude in the pupil after FFT - was complex64
        self.tempImgMem=numpy.zeros((nfft,nfft),numpy.float32)#self.fpDataType)# Square modulus of the FFT of pupilAmplitude - was float64
        self.binimgMem=None
        if self.nimgNotNfft:
            self.binimgMem=numpy.zeros((nimg,nimg),numpy.float32)#self.fpDataType)#binned image
        ##Allocation of the memory for the short exposure PSF
        ##Array storing the input phase
        if type(fpgaarr)==type(None):
            self.instImgMem=numpy.zeros((nimg,nimg),numpy.float32)#self.fpDataType)# Instantaneous image - dummy size since copied later - was float64
            self.phsMem=numpy.zeros((npup*npup*pupmult,),numpy.float32)#self.fpDataType)#was float64
            
        else:
            self.phsMem=fpgaarr
            self.instImgMem=fpgaarr
        #self.longExpPSFMem=Numeric.zeros((nfft,nfft),Numeric.Float64)# Long exposure image array
        self.fftTmpMem=numpy.zeros((nfft**2+10*nfft,),numpy.complex64)#self.cpDataType)#was complex64
        #now initialise the science objects with this memory.
        for this in self.thisObjList:
            sciObj=this.sciObj
            sciObj.initMem(self.fftTmpMem,self.pupilAmplitudeMem,self.focusAmplitudeMem,self.tempImgMem,self.instImgMem,self.phsMem,self.binimgMem)
            sciObj.initialiseFPGA(fpid=fpid,fpgaInfo=fpgaInfo)
            sciObj.initProfiles(self.control["useFPGA"])
            f=int(0.4*this.sciObj.nimg)
            t=int(0.6*this.sciObj.nimg)
            this.sciObj.instImgView=this.sciObj.instImg[f:t,f:t]
            this.sciObj.longExpPSFView=this.sciObj.longExpPSF[f:t,f:t]



            

            
    def prepareNextIter(self):
        
        this=self.thisObjList[self.currentIdObjCnt]
        for attr in dir(this):
            if attr not in ["__doc__","__module__","__init__","idstr"]:
                val=getattr(this,attr)
                setattr(self,attr,val)
    def endThisIter(self):
        """Copy things back to the this object..."""
        pass
##         this=self.thisObjList[self.currentIdObjCnt]
##         this.n_integn=self.n_integn
##         this.isamp=self.isamp
    def strParams(self,indx,iter=None):
        if iter==None:
            txt="{"
        else:
            txt="%d: {"%iter
        d=self.thisObjList[indx].sciObj.dictScience
        for k in d.keys():
            txt+="'%s': "%k
            try:
                txt+="%g, "%d[k]
            except:
                txt+="%s, "%d[k]
        txt+="}"#, RMS: %g"%self.thisObjList[indx].sciObj.phaseRMS
        return txt
        

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
    def __del__(self):
        self.endSim()

    def endSim(self):
        #print "Destroying science module:  Release of FFTW plan and FFTW associated data"
        #nrfft.fftw_destroy_plan(self.fftPlan)
        #nrfft.fftw_cleanup()
        if self.ended==0:
            self.ended=1
            idtxt=""
            if self.config.this.simID!="":
                idtxt=" (%s)"%self.config.this.simID
            print "INFORMATION: Science results for %s%s:"%(self.moduleName,idtxt)
            timestamp=time.strftime("%y%m%d_%H%M%S")
            for this in self.thisObjList:
                # compute the OTF, if not already done.
                if this.sciObj.computeOTF==0:
                    sl=numpy.fft.fftshift(this.sciObj.longExpPSF)
                    this.sciObj.dictScience['strehlOTF']=numpy.fft.fft2(sl,s=(sl.shape[0]*2,sl.shape[1]*2)).sum()/this.sciObj.diffOTFSum

                print "**%s**:"%(str(this.objID))
                maxLen=0
                for k in this.sciObj.dictScience:
                    if len(k)>maxLen: maxLen=len(k)
                for k in this.sciObj.dictScience:
                    print "  %s**%s** =%s"%(
                        " "*(maxLen-len(k)),k,str(this.sciObj.dictScience[k]) )
                rmstxt=""
                if this.sciObj.phaseRMScnt>0:
                    m=this.sciObj.phaseRMSsum/this.sciObj.phaseRMScnt
                    rmstxt=", RMS: %g +- %g +- %g"%(m,numpy.sqrt(this.sciObj.phaseRMSsum2/this.sciObj.phaseRMScnt-m*m),numpy.sqrt((this.sciObj.phaseRMSsum2/this.sciObj.phaseRMScnt-m*m)/this.sciObj.phaseRMScnt))#note, when plotting, plot rms with error bars of the smaller one - this is the error of the mean, not the sample distribution...
                if this.sciObj.fitsFilename!=None:
                    #if os.path.exists(self.fitsFilename):
                    #os.remove(self.fitsFilename);
                    self.mkdirForFile(this.sciObj.fitsFilename)
                    head=[]
                    head.append("NAME    = '%s'"%this.objID)
                    head.append("NINTEG  = %d"%this.sciObj.n_integn)
                    head.append("SCALE   = %g"%this.sciObj.pix_scale)
                    head.append("SAVESTR = '%s'"%this.sciObj.saveFileString)
                    head.append("SIMID   = '%s'"%self.config.this.simID)
                    head.append("SCISAMP = %d"%this.sciObj.sciPSFSamp)
                    head.append("BATCHNO = %d"%self.config.this.batchNumber)
                    head.append("TIMSTAMP= '%s'"%timestamp)
                    for key in this.sciObj.dictScience.keys():
                        head.append("%-8s= %s"%(key[:8],str(this.sciObj.dictScience[key])))
                    if rmstxt!="":
                        head.append("RMS     = %s"%rmstxt)
                    if len(this.sciObj.saveFileString)>0:
                        head.append("SIMINFO = '%s'"%this.sciObj.saveFileString)
                    util.FITS.Write(this.sciObj.longExpImg/this.sciObj.n_integn,this.sciObj.fitsFilename,extraHeader=head,writeMode="a",splitExtraHeader=1)
                if this.sciObj.sciFilename!=None:
                    self.mkdirForFile(this.sciObj.sciFilename)
                    f=open(this.sciObj.sciFilename,"a")
                    f.write("%s%s%s (%dx%d iters, batchno %d %s): %s%s\n"%(str(this.objID),this.sciObj.saveFileString,idtxt,this.sciObj.n_integn,this.sciObj.sciPSFSamp,self.config.this.batchNumber,timestamp,str(this.sciObj.dictScience),rmstxt))
                    f.close()
                if this.sciObj.histFilename!=None:
                    self.mkdirForFile(this.sciObj.histFilename)
                    head=[]
                    head.append("NAME    = '%s'"%this.objID)
                    head.append("NINTEG  = %d"%this.sciObj.n_integn)
                    head.append("SCALE   = %g"%this.sciObj.pix_scale)
                    head.append("SAVESTR = '%s'"%this.sciObj.saveFileString)
                    head.append("SIMID   = '%s'"%self.config.this.simID)
                    head.append("SCISAMP = %d"%this.sciObj.sciPSFSamp)
                    head.append("BATCHNO = %d"%self.config.this.batchNumber)
                    head.append("TIMSTAMP= '%s'"%timestamp)
                    k=str(this.sciObj.dictScience.keys())
                    k=k.replace("'",'').replace("[","").replace("]","")

                    head.append("KEYS    = '%s'"%k)
                    util.FITS.Write(this.sciObj.history[:,:this.sciObj.historyCnt],this.sciObj.histFilename,extraHeader=head,writeMode="a",splitExtraHeader=1)
                if this.sciObj.luckyFilename!=None:
                    self.mkdirForFile(this.sciObj.luckyFilename)
                    #Now write the history.
                    head=[]
                    head.append("NAME    = '%s'"%this.objID)
                    head.append("NINTEG  = %d"%this.sciObj.n_integn)
                    head.append("SCALE   = %g"%this.sciObj.pix_scale)
                    head.append("SAVESTR = '%s'"%this.sciObj.saveFileString)
                    head.append("SIMID   = '%s'"%self.config.this.simID)
                    head.append("SCISAMP = %d"%this.sciObj.sciPSFSamp)
                    head.append("BATCHNO = %d"%self.config.this.batchNumber)
                    head.append("LUCKSAMP= %d"%this.sciObj.luckyNSampFrames)
                    head.append("TIMSTAMP= '%s'"%timestamp)
                    k=str(this.sciObj.luckyHistoryKeys)
                    k=k.replace("'",'').replace("[","").replace("]","")

                    head.append("KEYS    = '%s'"%k)
                    util.FITS.Write(this.sciObj.luckyHistory[:,:this.sciObj.luckyHistoryCnt],this.sciObj.luckyFilename,extraHeader=head,writeMode="a",splitExtraHeader=1)
                if this.sciObj.luckyFile!=None:
                    self.mkdirForFile(this.sciObj.luckyImgFilename)
                    # we have been writing lucky images, so file is already open - need to finalise it (write correct header size)
                    #Actually we need to open a FITS file, write the header, copy data from luckyFile, and then finalise it, and close luckyFile (which is a tempfile and so gets deleted).
                    pos=this.sciObj.luckyFile.tell()
                    dtype=this.sciObj.luckyImg.dtype.char
                    shape=[pos/this.sciObj.luckyImg.itemsize/this.sciObj.luckyImgSize**2,this.sciObj.luckyImgSize,this.sciObj.luckyImgSize]
                    head=[]
                    head.append("NAME    = '%s'"%this.objID)
                    head.append("NINTEG  = %d"%this.sciObj.n_integn)
                    head.append("SCALE   = %g"%this.sciObj.pix_scale)
                    head.append("SAVESTR = '%s'"%this.sciObj.saveFileString)
                    head.append("SIMID   = '%s'"%self.config.this.simID)
                    head.append("SCISAMP = %d"%this.sciObj.sciPSFSamp)
                    head.append("BATCHNO = %d"%self.config.this.batchNumber)
                    head.append("LUCKSAMP= %d"%this.sciObj.luckyNSampFrames)
                    head.append("TIMSTAMP= '%s'"%timestamp)
                    txt=util.FITS.MakeHeader(shape,dtype,extraHeader=head,doByteSwap=this.sciObj.luckyByteswap,extension=os.path.exists(this.sciObj.luckyImgFilename),splitExtraHeader=1)
                    f=open(this.sciObj.luckyImgFilename,"a")
                    f.write(txt)
                    this.sciObj.luckyFile.seek(0)
                    for i in range(pos//8192):
                        f.write(this.sciObj.luckyFile.read(8192))
                    if pos%8192!=0:
                        f.write(this.sciObj.luckyFile.read(pos%8192))
                    this.sciObj.luckyFile.close()
                    this.sciObj.luckyFile=None
                    pos=2880-f.tell()%2880
                    if pos<2880:
                        f.write(" "*pos)
                    f.close()
                    

    def mkdirForFile(self,fname):
        d,f=os.path.split(fname)
        if len(d)>0:
            if not os.path.exists(d):
                os.makedirs(d)
##     def initRadialProfile(self):
##         """Creates and initialises arrays to fastly compute the radial profile and the encircled energy profile of the PSF
##         Must be first called before calling computeRadialProfileAndEncircledEnergy
##         """	
##         ####  We first create a map of the square of the distances to the center of the PSF image
##         ####  The map of distances is computed with the dist function (see aosim/utils/dist.py for help)
##         ####  We use the square because we have then integer indexes
##         r2=dist(self.nfft)
##         r2*=r2 ## square : we have integer numbers
##         self.nnRad=Numeric.argsort(r2.flat) ##we flatten the grid of distances and sort the values
##         r2r=Numeric.take(r2.flat,self.nnRad) ## we sort the grid of distances
##         self.xRad=Numeric.nonzero(dif_(r2r)) ##we look when the grid of distances change of value
##         self.tabNbPointsRad=dif_(self.xRad) ## number of points per bin
##         self.rRad=Numeric.take(Numeric.sqrt(r2r),self.xRad) ##radius (in pixels) giving the horizontal axis for radial and encircled energy profiles

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
##         profil=tabEnergyPerPixBin/self.tabNbPointsRad

##         #### Allocation of the return result
##         result=Numeric.zeros((2,len(tabEncircledEnergy)),typecode=psf.typecode())

##         #### We first store the radial profile in line 1 of the returned array
##         result[0,0]=psfSort[0]
##         result[0,1:,]=profil[:,] ##radial profile of PSF

##         #### We  store the encircled energy profile in line 2 of the returned array
##         result[1,:,]=tabEncircledEnergy[:,]
##         return result

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
    

##     def computeEnsquaredEnergy(self,psf):
##         """Computes ensquared energy profile
##         """
##         #### We flatten the diffraction limited PSF and sort the values
##         psfSort=Numeric.take(psf.flat,self.nnSquare)

##         #### We compute the encircled energy of the diffraction limited PSF 
##         tabEnsquaredEnergy=Numeric.take(Numeric.cumsum(psfSort),self.xSquare)

##         #### We convert the profile into psf type and return the array
##         result=tabEnsquaredEnergy.astype(psf.typecode())
##         return result


##     def makeData(self):
##         """Prepare the data for sending over socket.  ie turn it into an array.
##         With the histories, it would actally make more sense to have them always
##         stored in an array not a list, and when the array is full, remove history logarithmically
##         so that you always have info about the start...
##         This can be called by the GUI every cycle or whatever, returning the data...
##         """
##         self.plot_data['strehl']=Numeric.array(self.strehl_hist)# Send plot data on demand
##         self.plot_data['d50']=Numeric.array(self.d50_hist)
##         self.plot_data['fwhm']=Numeric.array(self.fwhm_hist)
##         #self.plot_data['inhex']=Numeric.array(self.inhex_hist) ##for hexagonal geometry : not used
##         self.plot_data['inbox']=Numeric.array(self.inbox_hist)
##         ## The piston/tip/tilt removal code has already been done in atmos
##         ## Indeed no, it has been commented (I (FA) WOULD LIKE TO KNOW WHY - PLEASE PUT COMMENTS INTO YOUR SOURCE WHEN YOU DO MODIFICATIONS !!!!!!!!!!!!!!!!!!!!)
##         ## So we compute the piston of the phase
##         ##pist=Numeric.sum(Numeric.sum(self.inputData*self.pupil.fn))/Numeric.sum(Numeric.sum(self.pupil.fn))
##         self.plot_data['sciphs']= (self.inputData)*self.pup ##phase seen by the science camera
##         self.plot_data['scilongexpimg']=self.longexp_img[nfft/2-nfft/8:nfft/2+nfft/8,
##                                                          nfft/2-nfft/8:nfft/2+nfft/8]
##         self.plot_data['scishortexpimg']= self.inst_img[nfft/2-nfft/4:nfft/2+nfft/4,
##                                                         nfft/2-nfft/4:nfft/2+nfft/4]

##     def needData(self,msg=None):
##         """Always need data to keep it running but will call getData
##         (maybe more than once) during generateNext if wfs integration time
##         is longer than timestep."""
##         return False

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

    def generateNext(self,msg=None):
        """science main loop... gets phase data and calculates psf
        etc.  Then stores these in a history.  If data is required for
        plots, should be obtained using the gui...  Each call to
        generateNext calls el_dm.next() etc the same number of times
        as wfs does (when doing wfs integration).  This ensures the
        simulation runs in sync.  If here, the integration length is
        different from that of the wfs, then have to bear that in
        mind."""
        ##print "Executing science.next method"
        ##we first copy the phase in input data to the self.inputData array, taking only the values into the pupil
        self.ended=0
        t1=time.time()
        if self.debug!=None:
            print "science: in generateNext() (debug=%s)"%str(self.debug)
        if self.generate==1:
            if self.newDataWaiting:
                if self.parent.dataValid==1:
                    self.inputData=self.parent.outputData
                    self.dataValid=1
                else:
                    print "INFORMATION: Science: waiting for data from DM, but not valid"
                    self.dataValid=0
            if self.dataValid:
                self.sciObj.doScienceCalc(self.inputData,self.control,self.thisiter)
        else:
            self.dataValid=0
        if self.debug!=None:
            print "Science: done generateNext() (debug=%s)"%str(self.debug)
        self.generateNextTime=time.time()-t1
        if self.currentIdObjCnt==len(self.idstr)-1:
            self.thisiter+=1
            if self.control["zero_science"]>0:
                self.control["zero_science"]-=1

##     def doScienceCalc(self):
##         ##we compute the piston and remove it into the pupil
##         if self.atmosPhaseType=="phaseonly":
##             pist=Numeric.sum(Numeric.sum(self.inputData*self.pup))/self.pupsum ##computation of the piston
##             Numeric.put(self.phs.flat,self.idxPup,Numeric.take(self.inputData.flat,self.idxPup)-pist) ##we remove the piston only from stuff in the pupil.
##         else:
##             raise Exception("science: todo: don't know how to remove piston")


##         nfft=self.nfft
##         t1=time.time()
##         if self.control["calcRMSFile"]!=None:
##             rms=self.calcRMS(self.phs,self.pup)
##             f=open(self.control["calcRMSFile"],"a")
##             f.write("%g\t%g\n"%(self.thisiter*self.tstep,rms))
##             f.close()

##         if self.control["science_integrate"]:# Integrate if shutter is open
##             if self.control["zero_science"]:# Zero science image
##                 self.longexp_img*=0.
##                 self.n_integn=0
##             else: ## integration of the PSF
##                 self.computeShortExposurePSF(self.phs) ##compute short exposure PSF
##                 self.longexp_img=self.longexp_img+self.inst_img	# Integrate long exposure
##                 self.n_integn+=1 ##We increment the number of integrations used to compute long exposure PSF
##             ##we increment the isamp counter
##             self.isamp+=1
##             if (self.isamp>=self.scinSamp): ## we compute and store scientific parameters and store the PSF
##                 ##we reset isamp to 0
##                 self.isamp=0
##                 ##if self.control["zero_science"]=0 : all the parameters in the dictScience dictionary are equal to 0 ; otherwise we compute them in the computeScientificParameters function
##                 if (self.control["zero_science"]==1):
##                     self.control["zero_science"]=0
##                     self.dictScience['FWHM']=0.
##                     self.dictScience['d50']=0.
##                     self.dictScience['inbox']=0.
##                     self.dictScience['strehl']=0.
##                 else: ## we call the function to compute science parameters as the PSF has not been zeroed
##                     self.computeScientificParameters()

##                 ## we update the history lists
##                 ## first the strehl
##                 if (len(self.strehl_hist)==self.scienceListsSize): ## we arrive at the last element of the science parameter  list : we remove the 1st element to add the new one at the end so that the list has always scienceListsSize elements
##                     vout=self.strehl_hist.pop(0)
##                 self.strehl_hist.append(self.dictScience['strehl'])

##                 ##then the fwhm
##                 if (len(self.fwhm_hist)==self.scienceListsSize): ## we arrive at the last element of the science parameter  list : we remove the 1st element to add the new one at the end so that the list has always scienceListsSize elements
##                     vout=self.fwhm_hist.pop(0)
##                 self.fwhm_hist.append(self.dictScience['FWHM'])

##                 ##then the d50
##                 if (len(self.d50_hist)==self.scienceListsSize): ## we arrive at the last element of the science parameter  list : we remove the 1st element to add the new one at the end so that the list has always scienceListsSize elements
##                     vout=self.d50_hist.pop(0)
##                 self.d50_hist.append(self.dictScience['d50'])

##                 ##then the ensquared energy in a 0.2 arcsec aperture
##                 if (len(self.inbox_hist)==self.scienceListsSize): ## we arrive at the last element of the science parameter  list : we remove the 1st element to add the new one at the end so that the list has always scienceListsSize elements
##                     vout=self.inbox_hist.pop(0)
##                 self.inbox_hist.append(self.dictScience['inbox'])

##                 if self.fitsFilename!=None:
##                     #self.psfWrittenToDisk[1]=self.longexp_img/self.n_integn ##we store the AO-corrected PSF
##                     if os.path.exists(self.fitsFilename): ##if the file already exists : we erase it
##                         os.remove(self.fitsFilename);
##                     util.FITS.Write(self.longexp_img/self.n_integn,self.fitsFilename)


##                 if self.write_to_file: ##do we want to write the arrays on disk
##                     f=open(self.filename,'w')# Science parameter file
##                     psf_params={}       # Save PSF data to file
##                     psf_params['scale']=self.pix_scale
##                     psf_params['strehl']=self.dictScience['strehl']
##                     ##psf_params['inhex0.3']=self.inhex
##                     psf_params['fwhm']=self.dictScience['FWHM']
##                     psf_params['d50']=self.dictScience['d50']
##                     ##psf_params['encirc']=self.encirc
##                     ##psf_params['azav']=self.azpsf
##                     ##psf_params['psf']=self.plot_data['scilongexpimg']
##                     cPickle.dump(psf_params,f)
##                     f.close()

##         ## For debugging purpose
##         if(self.timing):print thread.getIdent(),"science",time.time()-t1

##         #self.anyplots=0
##         #for plot in self.plots:
##         #  if(self.controlDict[plot]):# Is this plot requested ?
##         #    self.anyplots=1
##         #    for obj in self.objs:
##         #     if( int(obj) <= len(self.plot_source[plot]) ):    
##         #      if(self.controlDict[obj]):
##         #	if(rankDict[self.plot_source[plot][int(obj)-1]] == comm.rank):# Is this process the data source ?
##         #	  tag=100+int(self.plots.index(plot))# Assign a unique message tag
##         #	  comm.send(self.plot_data[plot],rankDict['plotter'],tag)
##         #print 'Science iteration done'
##         if self.debug!=None:
##             print "science: generateNext done (debug=%s)"%str(self.debug)
##         return None


## #### Computes the diffraction limited PSF
##     def computeDiffPSF(self):
##         """computeDiffPSF function : computes the diffraction limited PSF
##         Function written  by FA"""
##         if self.atmosPhaseType=="phaseamp":
##             self.computeShortExposurePSF([1,1])
##         else:
##             self.computeShortExposurePSF(1)
##         if self.keepDiffPsf:
##             self.diffPSF=self.inst_img.copy()
##         else:
##             self.diffPSF=self.inst_img#will be overwritten next time computeShortExposurePSF is called...
##         return self.diffPSF

##         ##We copy the arrays to have contigous arrays and convert it into Float64
##         #pup2=self.pup.copy()
##         #pup2=pup2.astype("d")
## ##         npup=self.npup

## ##         ##We fill the complex amplitude
## ##         self.pupilAmplitude[npup:,]=0.#clear the array (may get changed by fft)
## ##         self.pupilAmplitude[:npup,npup:,]=0.
## ##         self.pupilAmplitude.real[:npup,:npup]=self.pup ##*cos(phi)=1
## ##         self.pupilAmplitude.imag[:npup,:npup]=0. ##*sin(phi)=0
## ##         ##We call the FFTW function
## ##         nrfft.fftw_execute(self.fftPlan)
## ##         ##We compute the intensity in the focal plane
## ##         self.tempImg[:,]=(Numeric.absolute(self.focusAmplitude))**2
## ##         ##We compute the PSF by using fliparray2
## ##         self.diffPSF=fliparray2(self.tempImg)# Flip quadrants with definition consistent with FFT coordinate definition
## ##         self.diffPSF/=Numeric.sum(Numeric.sum(self.diffPSF)) ##We normalise the instantaneous PSF to 1


## #### Computes the short exposure PSF
##     def computeShortExposurePSF(self,phs):
##         """computeShortExposurePSF function : computes short exposure AO corrected PSF
##         Modifications made by FA"""

## ##        ##We copy the arrays to have contigous arrays and convert it into Float64
## ##        phs2=phs.copy()
## ##        phs2=phs2.astype("d")
## ##        self.phs2=phs2
## ##        pup2=self.pup.copy()
## ##        pup2=pup2.astype("d")
## ##        self.pup2=pup2
##         npup=self.npup ##to have the dimensions of the input phase array

##         ##We fill the complex amplitude
##         if self.atmosPhaseType=="phaseonly":
##             self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
##             self.pupilAmplitude[:npup,npup:,]=0.
##             self.pupilAmplitude.real[:npup,:npup]=self.pup*Numeric.cos(phs)
##             self.pupilAmplitude.imag[:npup,:npup]=self.pup*Numeric.sin(phs)
##         elif self.atmosPhaseType=="phaseamp":#phs[1] is amplitude, phs[0] is phase
##             self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
##             self.pupilAmplitude[:npup,npup:,]=0.
##             self.pupilAmplitude.real[:npup,:npup]=self.pup*Numeric.cos(phs[0])*phs[1]
##             self.pupilAmplitude.imag[:npup,:npup]=self.pup*Numeric.sin(phs[0])*phs[1]
##         elif self.atmosPhaseType=="realimag":#phs in real/imag already.
##             self.pupilAmplitude[npup:,]=0.#clear array (may get changed by fft)
##             self.pupilAmplitude[:npup,npup:,]=0.
##             self.pupilAmplitude[:npup,:npup]=self.pup*phs

            
##         ##We call the FFTW function
##         nrfft.fftw_execute(self.fftPlan)
##         ##We compute the intensity in the focal plane
##         self.tempImg[:,]=(Numeric.absolute(self.focusAmplitude))**2
##         ##print phs2.shape,phs2.typecode(),self.tempImg.shape,self.tempImg.typecode(),pup2.shape,pup2.typecode()
##         ##We compute the PSF by using fliparray2
##         #self.inst_img=fliparray2(self.tempImg)# Flip quadrants with definition consistent with FFT coordinate definition
##         fliparray2(self.tempImg,self.inst_img)# Flip quadrants with definition consistent with FFT coordinate definition
##         self.inst_img/=Numeric.sum(Numeric.sum(self.inst_img)) ##We normalise the instantaneous PSF to 1

        
## ### Computes scientific parameters on the long-exposureAO-corrected PSF #############
##     def computeScientificParameters(self):
##         """computeScientificParameters function : computes FWHM, 
##            Modifications made by FA"""
##         ## We do the average of the PSF and normalise it to 1
##         self.longExpPSF[:,]=self.longexp_img#/self.n_integn
##         self.longExpPSF/=Numeric.sum(Numeric.sum(self.longExpPSF)) ##We normalise the instantaneous PSF to 1
##         #window(2,wait=1)
##         #fma()
##         #pli(self.longExpPSF)
                
##         ## We calculate PSF parameters ###
##         ## We start by the Strehl Ratio
##         nfft=self.nfft
##         strehl=self.longExpPSF[nfft/2,nfft/2]/self.diffn_core_en
##         self.dictScience['strehl']=strehl
##         ##print "Strehl=%g"%(strehl)

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
##             print "ERROR in science.py - no non-zero elements... Oh well! (%s)"%self.objID
##             i1=0
##         x1=x[i1]
##         t=Numeric.less(y,y[0]/2)
##         try:
##             i2=Numeric.nonzero(t)[0]
##         except:
##             print "ERROR2 in science.py - no non-zero elements... Oh well! (%s)"%self.objID
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
##         apSize=0.2
##         t=Numeric.greater(x,apSize)
##         try:
##             i1=Numeric.nonzero(t)[0]
##         except:
##             print "ERROR5 in science.py - no non-zero elements... Oh well!"
##             i1=0

##         x1=x[i1]
##         t=Numeric.less(x,apSize)
##         try:
##             i2=Numeric.nonzero(t)[-1]
##         except:
##             print "ERROR6 in science.py - no non-zero elements... Oh well!"
##             i1=0
            
##         x2=x[i2]
##         if x2==x1:
##             x2+=1
##         ##we compute the equation of the line between those two points
##         a=(y[i2]-y[i1])/(x2-x1)
##         b=y[i2]-a*x2
##         ##we then compute the energy
##         inbox=a*apSize+b
##         ##print "E(0.2 arcsec)=%g"%(inbox)
##         self.dictScience['inbox']=inbox

#### Convolution by an instrumental PSF. Requires to set img
##     def conv(self,img):
##         """Convolve image with a PSF"""
##         temp=FFT.real_fft2d(img)
##         convimg=FFT.inverse_real_fft2d(temp*self.sci_psf_fft)
##         return convimg
    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        txt=""
        #if self.sentPlots==0:
        #for i in xrange(len(self.thisObjList)):
        this=self.thisObjList[self.sentPlotsCnt]
        if this.idstr==None or this.idstr=="":
            id=""
        else:
            id=" (%s)"%this.idstr
        txt+="""<plot title="View central psf%s" when="cmd" ret="feedback" texttype="1" wintype="mainwindow">\n<cmd>\nfeedback=1-%s.control['viewCentral']\n%s.control['viewCentral']=feedback\nso=%s.thisObjList[%d].sciObj\nf=int(so.nimg*0.4)\nt=int(so.nimg*0.6)\nif feedback:\n so.instImgView=so.instImg[f:t,f:t]\n so.longExpPSFView=so.longExpPSFView[f:t,f:t]\nelse:\n so.instImgView=so.instImg\n so.longExpPSFView=so.longExpPSF\n</cmd>\nbutton=feedback\n</plot>\n"""%(id,objname,objname,objname,self.sentPlotsCnt)
        txt+="""<plot title="Instantaneous image%s" cmd="data=%s.thisObjList[%d].sciObj.instImgView" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt)
        txt+="""<plot title="science phase%s" cmd="data=%s.thisObjList[%d].sciObj.phs" ret="data" type="pylab" dim="2" when="rpt" palette="gray" />\n"""%(id,objname,self.sentPlotsCnt)
        txt+="""<plot title="long PSF%s" cmd="data=%s.thisObjList[%d].sciObj.longExpPSFView" ret="data" type="pylab" dim="2" when="rpt" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt)
        txt+="""<plot title="Science params%s" cmd="data=%s.strParams(%d,ctrl.thisiter)" ret="data" when="rpt%d" texttype="1" wintype="ownwindow" textreplace="0"/>\n"""%(id,objname,self.sentPlotsCnt,this.sciObj.scinSamp*this.sciObj.sciPSFSamp)
        #txt+="""<plot title="Science params%s" cmd="data=str(%s.thisObjList[%d].sciObj.dictScience)+' RMS:'+str(%s.thisObjList[%d].sciObj.phaseRMS)" ret="data" when="rpt%d" texttype="1" wintype="ownwindow" textreplace="0"/>\n"""%(id,objname,self.sentPlotsCnt,objname,self.sentPlotsCnt,this.sciObj.scinSamp)
        txt+="""<plot title="Science history%s" cmd="data=%s.thisObjList[%s].sciObj.history[:,:%s.thisObjList[%s].sciObj.historyCnt]" ret="data" when="rpt%d" type="pylab" dim="2" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt,objname,self.sentPlotsCnt,this.sciObj.scinSamp*this.sciObj.sciPSFSamp)
        txt+="""<plot title="Zero science%s" cmd="%s.control['zero_science']+=1;data='Science zeroed'" ret="data" when="cmd" texttype="1" wintype="mainwindow"/>\n"""%(id,objname)
        txt+="""<plot title="Compute OTF Strehl%s" when="cmd" ret="feedback" texttype="1" wintype="mainwindow">\n<cmd>\nso=%s.thisObjList[%d].sciObj\nfeedback=1-so.computeOTF\nso.computeOTF=feedback\n</cmd>\nbutton=feedback\n</plot>\n"""%(id,objname,self.sentPlotsCnt)
        txt+="""<plot title="Lucky history%s" cmd="data=%s.thisObjList[%s].sciObj.luckyHistory[:,:%s.thisObjList[%s].sciObj.luckyHistoryCnt]" ret="data" when="rpt%d" type="pylab" dim="2" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt,objname,self.sentPlotsCnt,this.sciObj.sciPSFSamp*this.sciObj.luckyNSampFrames)
        txt+="""<plot title="Lucky image%s" cmd="data=%s.thisObjList[%s].sciObj.luckyLastImg" ret="data" when="rpt%d" type="pylab" dim="2" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt,this.sciObj.sciPSFSamp*this.sciObj.luckyNSampFrames)

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
        paramList.append(base.dataType.dataType(description="atmosPhaseType",typ="s",val="phaseonly",comment="What data do phase screens contain (ie amplitude too?)"))
        #paramList.append(base.dataType.dataType(description="sourceLam",typ="f",val="1650.",comment="Science source wavelength"))
        paramList.append(base.dataType.dataType(description="timing",typ="i",val="0",comment="extra timing info"))
        paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="scinfft",typ="eval",val="this.globals.npup*2",comment="FFT size for science image (oversample pupil)"))
        paramList.append(base.dataType.dataType(description="scinimg",typ="eval",val="this.globals.scinfft",comment="Image size for science image (oversample pupil, binned from FFT)"))
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(this.globals.npup,this.globals.ntel/2,this.globals.ntel/2*this.globals.telSec/this.globals.telDiam,TODO)",comment="telescope pupil"))
        paramList.append(base.dataType.dataType(description="dmpupil",typ="eval",val="this.globals.pupil.dmpupil",comment="DM pupil"))
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="8.",comment="Telescope diameter (m)"))
        #paramList.append(base.dataType.dataType(description="arcsecrad",typ="eval",val="Numeric.pi/180./3600.",comment="convert arc seconds to radians"))
        paramList.append(base.dataType.dataType(description="science_integrate",typ="i",val="1",comment="integrate the science image..."))
        paramList.append(base.dataType.dataType(description="zero_science",typ="i",val="0",comment="zero the science image"))
        #paramList.append(base.dataType.dataType(description="simFilename",typ="s",val="sci.cpickle",comment="name of file to which science results are saved."))
        paramList.append(base.dataType.dataType(description="scifitsFilename",typ="eval",val="None",comment="fits image filename"))
        paramList.append(base.dataType.dataType(description="scicsvFilename",typ="eval",val="None",comment="csv results filename"))
        paramList.append(base.dataType.dataType(description="scinSamp",typ="i",val="20",comment="Create science image parameters every this many iterations that a science image is created (ie this times sciPSFSamp)."))
        paramList.append(base.dataType.dataType(description="sciPSFSamp",typ="i",val="1",comment="Create science image every this many iterations."))
        paramList.append(base.dataType.dataType(description="hist_list_size",typ="i",val="100",comment="number of history elements to store"))
        return paramList


#### Adds an exception to stop the program and go to the prompt
if __name__=="__main__":
    batchno=0
    import base.readConfig
    config=base.readConfig.AOXml("params.xml",batchno=batchno)    
    s=science(None,config)
    raise Exception("Can't test like this")
    print "Doing nothing!"


