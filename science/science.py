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
        self.timing=config.getVal("timing",default=0) ##for debugging purpose
        self.doneFinalInit=0
        self.nimgNotNfft=0
        self.ended=1
        self.control={}
        self.sciOverview=config.getVal("sciOverview",raiseerror=0)
        if self.sciOverview==None:
            print("DEPRECATION: WARNING: Please use util.sci.sciOverview to describe science objects in config file")
            #leave control as global (shared between all resource sharing objts)
            self.control["zero_science"]=config.getVal("zero_science",default=10) ##idem - zero for 10 iterations...
            self.control["science_integrate"]=config.getVal("science_integrate",default=1)
            self.control["calcRMS"]=self.config.getVal("science_calcRMS",default=0)##idem
            self.control["viewCentral"]=0.4#if set, the gui will only be sent the central 20% of the psfs.
            self.control["lucky_integrate"]=config.getVal("lucky_integrate",default=0)
        self.thisiter=0
        self.sentPlotsCnt=0

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
        tstep=this.config.getVal("tstep")
        if self.sciOverview==None:
            apt=this.config.getVal("atmosPhaseType",default="phaseonly")
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
            usedmpup=this.config.getVal("usedmpup",default=0)
            if usedmpup:
                pup=numpy.array((pupil.fn*dmpupil),numpy.int8)#.astype(numpy.int8)# Pupil seen by Science cam
            else:
                pup=pupil.fn

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
                if luckyHistorySize==0:
                    luckyHistorySize=100
            luckyByteswap=this.config.getVal("luckyByteswap",default=0)
            saveFileString=this.config.getVal("scisaveFileString",raiseerror=0)#an optional string that can be used to identify a given simulation in the saved results... 
            diffPsfFilename=this.config.getVal("sciDiffPsfFilename",default=None,raiseerror=0)
            keepDiffPsf=this.config.getVal("keepDiffPsf",default=0)#overwrite mem?
            scienceListsSize=this.config.getVal("hist_list_size",default=this.config.getVal("AOExpTime")/tstep/scinSamp)
            if scienceListsSize==0:
                scienceListsSize=100
            inboxDiamList=this.config.getVal("inboxDiamList",default=[0.2])
            histFilename=this.config.getVal("histFilename",default=None,raiseerror=0)
            userFitsHeader=None
            psfEnergyToSave=0.
            psfMinSize=10
            nimgLongExp=None
        else:#use sciOverview.
            this.sciInfo=self.sciOverview.getSciByID(idstr)
            #a slight bodge - each will overwrite control...
            self.control["zero_science"]=this.sciInfo.zeroPsf
            self.control["science_integrate"]=this.sciInfo.integrate
            self.control["calcRMS"]=this.sciInfo.calcRMS
            self.control["viewCentral"]=0.4#if set, the gui will only be sent the central 20% of the psfs.
            luckyObj=this.sciInfo.luckyObj
            if luckyObj==None:
                self.control["lucky_integrate"]=0
            else:
                self.control["lucky_integrate"]=luckyObj.integrate
            apt=this.sciInfo.phaseType
            pupil=this.sciInfo.pupil ##telescope entrance pupil
            npup=pupil.shape[0]#this.config.getVal("npup") ##number of linear pixels for pupil
            nfft=this.sciInfo.nfft
            if nfft==None:
                nfft=npup*2#int(2**numpy.ceil(numpy.log2(2*npup))) ##fft size
            nimg=this.sciInfo.nimg
            if nimg==None:
                nimg=nfft
            if nimg!=nfft:
                self.nimgNotNfft=1
            realPupil=this.sciInfo.realPupil
            sciPath=this.sciInfo.sciPath
            if type(realPupil)!=type(None):
                realPupil=realPupil.fn
            try:
                dmpupil=pupil.dmpupil
            except:
                dmpupil=this.sciInfo.dmpupil
            usedmpup=this.sciInfo.usedmpup
            if usedmpup:
                pup=numpy.array((pupil.fn*dmpupil),numpy.int8)
            else:
                pup=pupil.fn

            fitsFilename=this.sciInfo.psfFilename
            sciFilename=this.sciInfo.summaryFilename
            sciPSFSamp=this.sciInfo.psfSamp#1
            scinSamp=this.sciInfo.nsamp#10
            saveFileString=this.sciInfo.saveString
            diffPsfFilename=this.sciInfo.diffPsfFilename
            keepDiffPsf=this.config.getVal("keepDiffPsf",default=0)#overwrite mem?
            scienceListsSize=this.sciInfo.histListSize
            if scienceListsSize==None:
                scienceListsSize=this.config.getVal("AOExpTime")/tstep/scinSamp
                if scienceListsSize==0:
                    scienceListsSize=100
            inboxDiamList=this.sciInfo.inboxDiamList
            userFitsHeader=this.sciInfo.userFitsHeader
            psfEnergyToSave=this.sciInfo.psfEnergyToSave
            psfMinSize=this.sciInfo.psfMinSize
            nimgLongExp=this.sciInfo.nimgLongExp
            histFilename=this.sciInfo.histFilename
            if luckyObj!=None:
                luckyFilename=luckyObj.filename
                luckyImgFilename=luckyObj.imgFilename
                luckyImgSize=luckyObj.imgSize
                luckyNSampFrames=luckyObj.nSampFrames
                luckyHistorySize=luckyObj.histSize
                if luckyHistorySize==None:
                    luckyHistorySize=int(this.config.getVal("AOExpTime")/(tstep*sciPSFSamp*luckyNSampFrames))
                    if luckyHistorySize==0:
                        luckyHistorySize=100
                luckyByteswap=luckyObj.byteswap
            else:
                luckyFilename=luckyImgFilename=luckyImgSize=None
                luckyNSampFrames=1
                luckyHistorySize=int(this.config.getVal("AOExpTime")/(tstep*sciPSFSamp*luckyNSampFrames))
                if luckyHistorySize==0:
                    luckyHistorySize=100
                luckyByteswap=0


        atmosGeom=this.config.getVal("atmosGeom",default=None,raiseerror=0)
        sci_lam=None
        phsLam=None
        if atmosGeom!=None:
            sci_lam=atmosGeom.sourceLambda(sciPath)
            phsLam=atmosGeom.phaseLambda(sciPath)
            telDiam=atmosGeom.telDiam
        else:#depreciated
            telDiam=this.config.getVal("telDiam") ##Telescope diameter 
        if sci_lam==None:
            sci_lam=this.config.getVal("sourceLam")##Imaging wavelength
            print "Warning - sci_lam not in atmosGeom - using %g"%sci_lam
        if phsLam==None:
            phsLam=sci_lam
        print "science %s - Using wavelength of %g nm (phase at %g nm)"%(idstr,sci_lam,phsLam)



        phaseMultiplier=phsLam/sci_lam
        asRad=numpy.pi/180./3600.#this.config.getVal("arcsecRad") ##Conversion rad-> arcsec
        L_D= ( (sci_lam*1.e-9)/telDiam )/asRad ##Diffraction limited resolution
        pix_scale=L_D*float(npup)/float(nfft)*float(nfft)/float(nimg) ##pixel scale (arcsec/pixel in the science image).  
        ##name of the FITS file to store the PSF
        #now create the science object that will do most of the work...
        this.sciObj=util.sci.science(npup,nfft,pup,nimg=nimg,atmosPhaseType=apt,tstep=tstep,keepDiffPsf=keepDiffPsf,pix_scale=pix_scale,fitsFilename=fitsFilename,diffPsfFilename=diffPsfFilename,scinSamp=scinSamp,sciPSFSamp=sciPSFSamp,scienceListsSize=scienceListsSize,debug=self.debug,timing=self.timing,allocateMem=0,realPup=realPupil,fpDataType=self.fpDataType,inboxDiamList=inboxDiamList,sciFilename=sciFilename,saveFileString=saveFileString,nthreads=self.nthreads,histFilename=histFilename,phaseMultiplier=phaseMultiplier,luckyNSampFrames=luckyNSampFrames,luckyFilename=luckyFilename,luckyImgFilename=luckyImgFilename,luckyImgSize=luckyImgSize,luckyHistorySize=luckyHistorySize,luckyByteswap=luckyByteswap,userFitsHeader=userFitsHeader,psfEnergyToSave=psfEnergyToSave,psfMinSize=psfMinSize,nimgLongExp=nimgLongExp)

        

        
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
        nimgmax=0
        nimgLongExpMax=0
        for this in self.thisObjList:
            if this.sciObj.nfft>nfftmax:
                nfftmax=this.sciObj.nfft
            if this.sciObj.npup>npupmax:
                npupmax=this.sciObj.npup
            if this.sciObj.atmosPhaseType in ["phaseamp","realimag"]:
                pupmult=2
            if this.sciObj.nimg>nimgmax:
                nimgmax=this.sciObj.nimg
            if this.sciObj.nimgLongExp>nimgLongExpMax:
                nimgLongExpMax=this.sciObj.nimgLongExp
        npup=npupmax
        nfft=nfftmax
        nimg=nimgmax
        nimgLongExp=nimgLongExpMax
        #now allocate memory at the max size needed.  These memory arrays will then be shared between all the sciObj objects.

        self.pupilAmplitudeMem=numpy.zeros((nfft,nfft),numpy.complex64)#self.cpDataType) ##Complex amplitude in the pupil - was complex64
        self.focusAmplitudeMem=numpy.zeros((nfft,nfft),numpy.complex64)#self.cpDataType) ##Complex amplitude in the pupil after FFT - was complex64
        self.tempImgMem=numpy.zeros((nfft,nfft),numpy.float32)#self.fpDataType)# Square modulus of the FFT of pupilAmplitude - was float64
        self.binimgMem=None
        if self.nimgNotNfft:
            self.binimgMem=numpy.zeros((nimg,nimg),numpy.float32)#self.fpDataType)#binned image
        ##Allocation of the memory for the short exposure PSF
        ##Array storing the input phase
        self.instImgMem=numpy.zeros((nimg,nimg),numpy.float32)#self.fpDataType)# Instantaneous image - dummy size since copied later - was float64
        self.phsMem=numpy.zeros((npup*npup*pupmult,),numpy.float32)#self.fpDataType)#was float64
        #self.longExpPSFMem=Numeric.zeros((nfft,nfft),Numeric.Float64)# Long exposure image array
        self.fftTmpMem=numpy.zeros((nfft**2+10*nfft,),numpy.complex64)#self.cpDataType)#was complex64
        self.longExpPSFMem=numpy.zeros((nimgLongExp,nimgLongExp),self.fpDataType)
        #now initialise the science objects with this memory.
        for this in self.thisObjList:
            sciObj=this.sciObj
            sciObj.initMem(self.fftTmpMem,self.pupilAmplitudeMem,self.focusAmplitudeMem,self.tempImgMem,self.instImgMem,self.phsMem,self.binimgMem,self.longExpPSFMem)
            sciObj.initProfiles()
            f=int(0.4*this.sciObj.nimg)
            t=int(0.6*this.sciObj.nimg)
            this.sciObj.instImgView=this.sciObj.instImg[f:t,f:t]
            #this.sciObj.longExpPSFView=this.sciObj.longExpPSF[f:t,f:t]



            

            
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
            try:#try to get the git commit label...
                commit=os.popen("(cd %s && git log -1)"%os.path.split(__file__)[0]).read().split("\n")[0][7:]
            except:
                commit="NOT FOUND"
                print "git commit info not found"
            
            
            for this in self.thisObjList:
                #final computation of science parameters:
                this.sciObj.longExpPSF[:]=this.sciObj.longExpImg/(numpy.sum(this.sciObj.longExpImg)+this.sciObj.clippedEnergy)##We normalise the instantaneous PSF to 1
                clippedFrac=this.sciObj.clippedEnergy/(this.sciObj.longExpImg.sum()+this.sciObj.clippedEnergy)
                this.sciObj.computeScientificParameters()

                # compute the OTF, if not already done.
                #if this.sciObj.computeOTF==0:
                #    sl=numpy.fft.fftshift(this.sciObj.longExpPSF)
                #    this.sciObj.dictScience['strehlOTF']=numpy.fft.fft2(sl,s=(sl.shape[0]*2,sl.shape[1]*2)).sum()/this.sciObj.diffOTFSum

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
                    print rmstxt
                if this.sciObj.fitsFilename!=None:
                    #if os.path.exists(self.fitsFilename):
                    #os.remove(self.fitsFilename);
                    self.mkdirForFile(this.sciObj.fitsFilename)
                    head=[]
                    head.append("NAME    = '%s'"%this.objID)
                    head.append("NINTEG  = %d"%this.sciObj.n_integn)
                    head.append("SCALE   = %g /pxl scale in arcsec/pxl"%this.sciObj.pix_scale)
                    head.append("SAVESTR = '%s'"%this.sciObj.saveFileString)
                    head.append("SIMID   = '%s'"%self.config.this.simID)
                    head.append("SCISAMP = %d"%this.sciObj.sciPSFSamp)
                    head.append("BATCHNO = %d"%self.config.this.batchNumber)
                    head.append("TIMSTAMP= '%s'"%timestamp)
                    head.append("COMMIT  = '%s'"%commit)
                    head.append("CLIPPEDE= %g"%clippedFrac)
                    for key in this.sciObj.dictScience.keys():
                        head.append("%-8s= %s"%(key[:8],str(this.sciObj.dictScience[key])))
                    if rmstxt!="":
                        head.append("RMS     = %s"%rmstxt)
                    if len(this.sciObj.saveFileString)>0:
                        head.append("SIMINFO = '%s'"%this.sciObj.saveFileString)
                    if this.sciObj.userFitsHeader is not None:
                        head+=this.sciObj.userFitsHeader
                    #img=this.sciObj.longExpImg/this.sciObj.n_integn
                    img=this.sciObj.longExpPSF
                    if this.sciObj.psfEnergyToSave!=0:#shrink the psf until it contains this fraction of energy.  e.g. 0.99
                        img,excludedEnergy=self.clipPsfToFluxFrac(img,this.sciObj.psfEnergyToSave,this.sciObj.psfMinSize)
                        head.append("EExluded= %g"%excludedEnergy)
                    util.FITS.Write(img,this.sciObj.fitsFilename,extraHeader=head,writeMode="a",splitExtraHeader=1)
                if this.sciObj.sciFilename!=None:
                    self.mkdirForFile(this.sciObj.sciFilename)
                    f=open(this.sciObj.sciFilename,"a")
                    f.write("%s%s%s (%dx%d iters, batchno %d %s commit %s): %s%s\n"%(str(this.objID),this.sciObj.saveFileString,idtxt,this.sciObj.n_integn,this.sciObj.sciPSFSamp,self.config.this.batchNumber,timestamp,commit,str(this.sciObj.dictScience),rmstxt))
                    f.close()
                if this.sciObj.histFilename!=None and this.sciObj.history!=None:
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
                    head.append("COMMIT  = '%s'"%commit)
                    k=str(this.sciObj.dictScience.keys())
                    k=k.replace("'",'').replace("[","").replace("]","")

                    head.append("KEYS    = '%s'"%k)
                    if this.sciObj.userFitsHeader is not None:
                        head+=this.sciObj.userFitsHeader
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
                    head.append("COMMIT  = '%s'"%commit)
                    k=str(this.sciObj.luckyHistoryKeys)
                    k=k.replace("'",'').replace("[","").replace("]","")

                    head.append("KEYS    = '%s'"%k)
                    if this.sciObj.userFitsHeader is not None:
                        head+=this.sciObj.userFitsHeader
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
                    head.append("COMMIT  = '%s'"%commit)
                    if this.sciObj.userFitsHeader is not None:
                        head+=this.sciObj.userFitsHeader
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
                    
    def clipPsfToFluxFrac(self,img,frac,minSize=10):
        size=minSize/2
        if size<1:
            size=1
        totEnergy=img.sum()
        energy=totEnergy*frac
        my=img.shape[0]//2
        mx=img.shape[1]//2
        tot=img[my-size:my+size,mx-size:mx+size].sum()
        while size<img.shape[0]//2 and tot<energy:
            size+=1
            tot+=img[my-size,mx-size:mx+size].sum()
            tot+=img[my+size-1,mx-size:mx+size].sum()
            tot+=img[my-size+1:my+size-1,mx-size].sum()
            tot+=img[my-size+1:my+size-1,mx+size-1].sum()
        img=img[my-size:my+size,mx-size:mx+size]
        excludedEnergy=totEnergy-img.sum()
        return img,excludedEnergy
                    
    def mkdirForFile(self,fname):
        d,f=os.path.split(fname)
        if len(d)>0:
            if not os.path.exists(d):
                os.makedirs(d)

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

    def getLongPsf(self,objno=0):
        sciObj=self.thisObjList[objno].sciObj
        img=sciObj.longExpImg/(sciObj.longExpImg.sum()+sciObj.clippedEnergy)#sciObj.n_integn
        if self.control["viewCentral"]!=0:
            f=int(sciObj.nimgLongExp*self.control["viewCentral"])
            img=img[f:-f,f:-f]
        return img


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
        #since the inst image isn't stored, what the plotter will see will be the last one computed.
        if self.thisObjList[-1].idstr==None or self.thisObjList[-1].idstr=="":
            id0=""
        else:
            id0=" (%s)"%self.thisObjList[-1].idstr
        txt+="""<plot title="Instantaneous image%s" cmd="data=%s.thisObjList[%d].sciObj.instImgView" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id0,objname,self.sentPlotsCnt)
        txt+="""<plot title="long PSF%s" cmd="data=%s.getLongPsf(%d)" ret="data" type="pylab" dim="2" when="rpt" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt)

        txt+="""<plot title="Science params%s" cmd="data=%s.strParams(%d,ctrl.thisiter)" ret="data" when="rpt%d" texttype="1" wintype="ownwindow" textreplace="0"/>\n"""%(id,objname,self.sentPlotsCnt,this.sciObj.scinSamp*this.sciObj.sciPSFSamp)
        txt+="""<plot title="View central psf%s" when="cmd" ret="feedback" texttype="1" wintype="mainwindow">\n<cmd>
if %s.control['viewCentral']==0:
 %s.control['viewCentral']=0.4
 feedback=1
else:
 %s.control['viewCentral']=0.
 feedback=0
so=%s.thisObjList[%d].sciObj
f=int(so.nimg*0.4)
t=int(so.nimg*0.6)
if feedback:
 so.instImgView=so.instImg[f:t,f:t]
else:
 so.instImgView=so.instImg
</cmd>\nbutton=feedback\n</plot>\n"""%(id,objname,objname,objname,objname,self.sentPlotsCnt)
        
        txt+="""<plot title="science phase%s" cmd="data=%s.thisObjList[%d].sciObj.phs" ret="data" type="pylab" dim="2" when="rpt" palette="gray" />\n"""%(id,objname,self.sentPlotsCnt)
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


