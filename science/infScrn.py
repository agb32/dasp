#todo: could be made more efficient by adding all new columns/rows at once,
#rather than shifting, adding, shifting, adding... this would save on shifts.
#Not sure how much of a speedup this would give.
#todo: think about whether its worth passing the whole phasescreen, rather than
#just the new rows each time - if the infAtmos object is on the same node, this
#will save computation in infAtmos, and take no time to pass.  However, if on
#a different node, will take longer.  Maybe, we need some option for this.
#Maybe outputData should contain both, which can then be selected.  Or maybe
#infAtmos should look at parent.screen to see if it exists, and if so, use it.
#Yes, thats probably the most elegant way...

##Infinite phase screen module
"""FAscrn.py : defines the FAscrn class implementing
Francois Assemat's infinite phase screens method
for phase screen computation

Version 1 : Adds columns in the horizontal direction, removes left ones
and add new ones on the right

Version 2 : Adds columns and rows in any direction

Version 3: Integration into AO-sim framework

cf Assemat et al., Optics Express, Vol. 14, Issue 3, pp. 988-999 (February 2006)
Things to think about:
Need 2 modules:
1 creates new rows/columns of atmospheric phase each iteration.  No
interpolation is done.  The number of columns/rows added depends on
wind velocity and direction.  If a non-integer number is required each
time, should do the following:
For example if 2.3 rows need to be added:
iter desired number  actual needed  number added
1    2.3             3 (ceil(2.3))  3
2    4.6             5 (ceil(4.6))  2
3    6.9             7 (ceil(6.9))  2
4    9.2             10    etc      3
5    11.5            12             2
In this case, you would always send 3 rows the the child (successor)
object.  However, sometimes, only 2 of these would contain valid data.


Other module is specific for a given target.
It takes the new rows/columns to get the translated (no
interpolation yet) phase (it should know when all or all-1 columns are
valid), selects the part of this that are relevent for the
star in question.  It then does interpolation to create the atmos
phase between this object and the pupil.  This can be sent to the DM/wfs.

During initial iterations, the small number of rows/columns is sent
each iteration until enough have been sent to create the whole phase.
This avoids having to send a different sized array on initialisation.

Note this there is a bug due to (I think - agb - ) numerical precision (float64)
that means that the resolution of the phasescreen shouldn't be too high.  eg don't try using a
1024x1024 phasescreen for a 8m telescope as it will blow up after adding 10 or so rows/columns.
If you go larger than this, then if may well not be able to invert parts of the phase covariance
matrix, and so will crash.  Be warned!

"""

import base.aobase,base.dataType,util.getNewCols#,cmod.utils

##for debug : plotting package
import time
import os
#import gist
import traceback
##Numeric Python import
import numpy
na=numpy

##random number generation function
import numpy.random as ra

##import required packages from scipy
##import gamma and modified bessel function
from scipy.special import kv,gamma

##import FFT functions (FFTW)
import scipy.fftpack as FFT

##import Linear Algebra functions (LAPACK)
import scipy.linalg as LA

##import the matrix vector function (BLAS)
#import util.matrix as matrix
import util.FITS
import cmod.scrn

## for fast MVM:
import util.dot as quick # will replace numpy.dot when I implement it for double

matrixdot=quick.dot#matrix.dot

def computeScrnSizeOld(thetas,phis,ntel,npup,telDiam,altitude):
    """The initial phase screen must be large enough to hold the phase for all objects viewed at a given time - ie including off axis objects.  However, even for an object 1 arcmin off axis at 10km, they are only 2.9m apart, so shouldn't be huge.
    thetas is a dictionary of source theta directions (arcseconds).
    phis is dictionary of source phi directions (degrees)
    altitude is dictionary of layer altitudes.
    This is used in the parameter file.
    """
    arcsecrad=2*na.pi/360./3600.
    degrad=2*na.pi/360.
    scrnSize={}
    for altkey in altitude.keys():
        xposlist=[]
        yposlist=[]
        layerAlt=altitude[altkey]
        for key in thetas.keys():
            xposlist.append(layerAlt*na.fabs(na.tan(thetas[key]*arcsecrad)*na.cos(phis[key]*degrad)))
            yposlist.append(layerAlt*na.fabs(na.tan(thetas[key]*arcsecrad)*na.sin(phis[key]*degrad)))
        maxx=max(xposlist)
        maxy=max(yposlist)
        scrnXPxls=int(na.ceil(npup+maxx*float(ntel)/telDiam))+1
        scrnYPxls=int(na.ceil(npup+maxy*float(ntel)/telDiam))+1
        scrnSize[altkey]=(scrnXPxls,scrnYPxls)
    print "scrnSize: %s"%str(scrnSize)
    return scrnSize
    
def computeScrnSize(thetas,phis,ntel,npup,telDiam,altitude,zenith):
    """Zenith in degrees"""
    zenith=0.
    arcsecrad=2*na.pi/360./3600.
    degrad=2*na.pi/360.
    scrnSize={}
    for altkey in altitude.keys():
        xposlist=[]
        yposlist=[]
        layerAlt=altitude[altkey]
        for key in thetas.keys():
            xposlist.append(layerAlt/numpy.cos(zenith*degrad)*na.fabs(na.tan(thetas[key]*arcsecrad)*na.cos(phis[key]*degrad)))
            yposlist.append(layerAlt*na.fabs(na.tan(thetas[key]*arcsecrad)*na.sin(phis[key]*degrad)))
        maxx=max(xposlist)
        minx=min(xposlist)
        maxy=max(yposlist)
        miny=min(yposlist)
        scrnXPxls=int(na.ceil(maxx*ntel/telDiam+npup+na.ceil(minx*ntel/telDiam))+1)
        scrnYPxls=int(na.ceil(maxy*ntel/telDiam+npup+na.ceil(miny*ntel/telDiam))+1)
        scrnSize[altkey]=(scrnXPxls,scrnYPxls)
    print "scrnSize: %s"%str(scrnSize)
    return scrnSize

def distanceMap(n,m=None,dy=0,dx=0,natype=na.float32):
    """
    The distanceMap function creates a rectangular array in which the value of
    each element is proportional to its frequency. This array may be used
    for a variety of purposes, including frequency-domain filtering and
    making pretty pictures.

    Syntax
        Result = dist(n,[m,dy,dx])


Arguments
    n  The number of lines in the resulting array.

    m  The number of columns in the resulting array. If M is omitted,
    the resulting array will be N by N.

    dy The offset from the offset in the vertical direction
    dx The offset from the offset in the horizontal direction
    dx The offset in y of the center
    naType : data type of the output array (Float32 by default)
    """
    #dimensions of the output array
    if m is None:
        m=n

    ## generation of the required axes             
    axe_x=(na.arange(m)-m/2.-dx)*1.0
    axe_y=(na.arange(n)-n/2.-dy)*1.0

    ## Creation of the grid of distances
    f=na.sqrt((axe_y[:,na.newaxis])**2.+(axe_x[na.newaxis,:])**2.)
    f=f.astype(natype)

    ##we return the grid
    return f
def computeInitialScreen(config,idstr=None):
    """computes the initial phase screen with FFT technique.
    This is not part of the class, so that it can be called by infAtmos also,
    to get exactly same starting point for phase... (since only the new rows
    are ever passed!).
    Note, the phase screen is actually 1 pxl bigger in each direction than dpix
    because of the way cols/rows will be added.  This extra row/col is all zero
    and is at either start or end, depending on wind direction.
    Also, not square - rectangular instead.

    """
    ##we first compute the physical size of the array required to perform the FFT
    so=config.searchOrder
    if idstr!=None:
        searchOrder=["infScrn_"+idstr,"infScrn","globals"]
    else:
        searchOrder=["infScrn","globals"]
    config.setSearchOrder(searchOrder)
    Dtel=config.getVal("telDiam")#diameter in m.
    dpix=config.getVal("npup")#diameter in pixels.
    atmosGeom=config.getVal("atmosGeom",default=None,raiseerror=0)
    if atmosGeom==None:
        scrnXPxls=config.getVal("scrnXPxls")#this must be large enough so that it can hold the phase for all objects viewed at a given time - ie including off axis objects.  However, even for an object 1 arcminute off axis, at 10km, they are viewed only 2.9m apart, so this won't be too huge... depreciated
        scrnYPxls=config.getVal("scrnYPxls")#depreciated
        seed=config.getVal("initPhsSeed")#depreciated
        strLayer=config.getVal("strength")#depreciated
        windDirection=config.getVal("windDirection")#in degrees depreciated
        vWind=config.getVal("vWind")#in m/s depreciated
        globR0=config.getVal("r0")
        L0=config.getVal("l0")
    else:
        scrnXPxls=atmosGeom.getScrnXPxls(idstr)
        scrnYPxls=atmosGeom.getScrnYPxls(idstr)
        seed=atmosGeom.layerInitSeed(idstr)
        strLayer=atmosGeom.layerStrength(idstr)
        windDirection=atmosGeom.layerWind(idstr)
        vWind=atmosGeom.layerSpeed(idstr)
        globR0=atmosGeom.r0
        L0=atmosGeom.l0
    print "Generating phasescreen %s size %s"%(str(idstr),str((scrnXPxls,scrnYPxls)))
    tstep=config.getVal("tstep")
    natype="d"#config.getVal("dataType")
    phaseArray=makeInitialScreen(dpix,Dtel,L0,scrnXPxls,scrnYPxls,seed,tstep,globR0,strLayer,natype,windDirection,vWind)
    ##we return the array corresponding to the pupil
    config.setSearchOrder(so)
    return phaseArray

def makeInitialScreen(dpix=1024,Dtel=42.,L0=30.,scrnXPxls=None,scrnYPxls=None,seed=0,tstep=0.05,globR0=0.2,strLayer=1.,natype=na.float64,windDirection=0.,vWind=10.):
    """dpix is the telescope aperture diameter, and dtel is tel diameter.
    The actual number of pixels used is scrnXPxls x scrnYPxls.
    """

    if scrnXPxls==None:
        scrnXPxls=dpix
    if scrnYPxls==None:
        scrnYPxls=dpix
    if seed!=None:
        try:
            fname="scrn/scrn%dD%gL%gx%dy%ds%dt%gr%gs%gd%gv%g.fits"%(dpix,Dtel,L0,scrnXPxls,scrnYPxls,seed,tstep,globR0,strLayer,windDirection,vWind)
        except:
            fname=None
            print "Could not make filename for screen - so generating screen..."
    else:#random seed
        fname=None
    if fname!=None and os.path.exists(fname):
        print "Loading existing screen: %s"%fname
        phaseArray=util.FITS.Read(fname)[1]
        return phaseArray

    scrnPxls=max(scrnXPxls,scrnYPxls)
    pixScale=Dtel/float(dpix)
    ro=globR0*(strLayer**(-3./5.))##we compute the ro in the considered layer
    colAdd=-vWind*na.cos(windDirection*na.pi/180)/pixScale*tstep#number of pixels to step each iteration (as float).
    rowAdd=-vWind*na.sin(windDirection*na.pi/180)/pixScale*tstep#number of pixels to step each iteration (as float).


    if seed==None:
        print "ERROR (possibly): computeInitialScreen - seed is None, so timer will be used, meaning that the initial screen cannot be replicated, so if both infScrn and infAtmos try to create, you will get a bad phasescreen.  If you wish to use a random seed, use int(time.time()) in the parameter file - though this will only work if all running in the same process."
    if L0>=Dtel:
        scrnSize=2*L0
    elif Dtel>=2*L0:
        scrnSize=Dtel
    else:
        scrnSize=2*L0
    ##size in pixels of the required phase screen
    nfft=int(na.around(dpix*1.*scrnSize/Dtel))#agb: scrnPxls was dpix
    nfft=max(nfft,scrnPxls)#now choose the maximum size...
    print "Initial phase screen size %d"%nfft
    ##gaussian standard 2D random variable
    ra.seed(seed)
    rn=ra.randn(nfft,nfft)

    ##creation of Von Karman spectrum with power spectrum in (f**2+fo**2)**(-11/6.)
    ##circular grid
    axe_x=na.arange(nfft)-nfft/2.
    f=na.sqrt((axe_x[:,na.newaxis])**2.+(axe_x[na.newaxis,:])**2.)

    #for FFT computations
    ff = FFT.fftshift(f) # to have 0 frequency in lower left corner
    fo=1./L0; #f0 = 1/L0

    ##factor in the power spectrum (0.023)
    fact=(gamma(11./6.))**2.;
    fact/=2*na.pi**(11./3.);
    fact*=((24./5)*gamma(6./5))**(5./6);

    ##f0 in pixels for power spectrum (cf FA's thesis, page 327 for explanation)
    phozero=fo*1.*scrnSize;

    ##we compute the phase power spectrum
    f = (ff**2.+phozero**2.)**(-11./12.)
    f*=na.sqrt(fact*((scrnSize/ro)**(5./3.))*(nfft*1.)**2.) ## Von Karman phase power spectrum : 0.023/ro*(f^2+fo^2)^(-11/6)

    ##FFT of white noise
    ff = FFT.fft2(rn);

    ## coloration of white noise by Von Karman spectrum
    ff.real*=f;
    ff.imag*=f

    ##inverse 2D Fourier Transform
    phi = (FFT.ifft2(ff)).real
    #phi=phi.astype(natype)

    #if colAdd is less than zero, we are adding now columns on the right of the array.
    #If rowAdd is less than zero, we are adding new rows at the bottom of the array (note, that if plotting in Gist, this is the top of the array).

    phaseArray=na.zeros((scrnYPxls+1,scrnXPxls+1),natype)#was dpix agbc swapped
    if colAdd<0:
        if rowAdd<0:
            phaseArray[1:,1:]=phi[:scrnYPxls,:scrnXPxls]#was dpix
        else:
            phaseArray[:-1,1:]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    else:
        if rowAdd<0:
            phaseArray[1:,:-1]=phi[:scrnYPxls,:scrnXPxls]#was dpix
        else:
            phaseArray[:-1,:-1]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    if fname!=None:
        if not os.path.exists("scrn/"):
            os.mkdir("scrn")
        util.FITS.Write(phaseArray,fname)
    return phaseArray


class infScrn(base.aobase.aobase):
    """
    Class infinitePhaseScreen : computes infinite phase screens
    using Francois Assemat's infinite phase screens method
    The created phase screen is rectuangular.
    """
	
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """config is an readConfig.AOXml object loaded with the correct
        configuration file"""
        ##extract data from XML configuration file
        base.aobase.aobase.__init__(self,None,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.degRad=na.pi/180.#self.config.getVal("degRad")
        self.niter=0
        self.nbColToAdd=1
        ##we first extract basic parameters
        ##Physical size of the telescope
        self.Dtel=config.getVal("telDiam")
        ##Number of pixels used to simulate the telescope pupil
        self.dpix=config.getVal("npup")
        # Total size of phasescreen (will be larger than dpix if there are off axis sources, and will also depend on height of the layer).
        self.atmosGeom=config.getVal("atmosGeom")
        self.scrnXPxls=self.atmosGeom.getScrnXPxls(self.idstr[0])
        self.scrnYPxls=self.atmosGeom.getScrnYPxls(self.idstr[0])
        self.windDirection=self.atmosGeom.layerWind(self.idstr[0])
        self.vWind=self.atmosGeom.layerSpeed(self.idstr[0])
        self.strLayer=self.atmosGeom.layerStrength(self.idstr[0])
        ##Physical size in meters per pixel
        self.pixScale=self.Dtel/float(self.dpix)
        ##Precision of the output screen array
        self.dataType="d"#config.getVal("dataType")
        if self.dataType not in ["d",na.float64]:
            print "WARNING - infScrn - it is not known whether statistics for infScrn will be right if not Float64.  Would be good if someone could test this... ie is it true Kolmogorov/von Karman"
        #self.altitude=config.getVal("layerAltitude")#in m.
        self.windDirRad=self.windDirection*self.degRad
        self.tstep=config.getVal("tstep")
        self.useCmodule=config.getVal("useCmodule",default=1)
        # if colAdd<zero, we are adding now columns on the right of the array.
        # If rowAdd<zero, we are adding new rows at the bottom of the array (note, that if plotting in Gist, this is the top of the array).
        self.colAdd=-self.vWind*na.cos(self.windDirRad)/self.pixScale*self.tstep#number of pixels to step each iteration (as float).
        self.rowAdd=-self.vWind*na.sin(self.windDirRad)/self.pixScale*self.tstep#number of pixels to step each iteration (as float).
        #This then tells us the size of the output array -
        self.maxColAdd=int(na.ceil(na.fabs(self.colAdd)))
        self.maxRowAdd=int(na.ceil(na.fabs(self.rowAdd)))
        self.sendWholeScreen=self.config.getVal("sendWholeScreen",default=0)
        if self.sendWholeScreen==0:
            # 2 arrays are returned, packed into 1.  The dimensions of these 2 arrays are: (maxColAdd, scrnYPxls) and (maxRowAdd, scrnXPxls)
            self.dims=(self.maxColAdd*(self.scrnYPxls+1)+self.maxRowAdd*(self.scrnXPxls+1),)#was dpix
        else:
            self.dims=(self.scrnYPxls+1,self.scrnXPxls+1)

        ##Test for the future config GUI
        if forGUISetup==1:
            ##specify characteristics of the output data for future "LabView like" GUI
            self.outputData=(self.dims,self.dataType)##specify output for GUI
        else: ##we continue the initialisation
            ##we extract the outer scale, the global Fried Parameter
            ##and the turbulence strength in the layer
            if self.atmosGeom==None:
                self.L0=config.getVal("l0")
                self.globR0=config.getVal("r0")
            else:
                self.L0=self.atmosGeom.l0
                self.globR0=self.atmosGeom.r0
            self.strLayerToPowMinusThreeOverFive=self.strLayer**(-3./5)
            self.stepFun=config.getVal("stepFunction",default=None,raiseerror=0)#function giving phase steps... as function of time.
            self.xstep=0.
            self.ystep=0.
            self.r0Fun=config.getVal("r0Function",default=None,raiseerror=0)#function giving r0 as function of time.
            ##we compute the ro in the considered layer
            if self.r0Fun!=None:
                self.computeR0()
            else:
                self.ro=self.globR0*(self.strLayerToPowMinusThreeOverFive)
            
            ##Number of previous rows or columns used to compute the new one (default value : 2)
            self.nbCol=2#config.getVal("nbCol",default=2)
            ##we extract the random number generator seed
            self.seed=config.getVal("seed",default=None)
            self.saveCovMat=config.getVal("saveInfPhaseCovMatrix",default=0)
            # Create the objects used to tell us how many extra cols/rows to add
            self.newCols=util.getNewCols.getNewCols(na.fabs(self.colAdd))
            self.newRows=util.getNewCols.getNewCols(na.fabs(self.rowAdd))

            ##we go now through the creation of the required matrices
            ##we compute first the phase covariance matrices
            if not os.path.exists("scrn/"):
                os.mkdir("scrn")
            fname="scrn/infScrnData%d_%g_%g_%d_%d.fits"%(self.scrnXPxls+1,self.L0,self.pixScale,self.nbColToAdd,self.nbCol)
            covMatPhix=None
            if os.path.exists(fname):
                print "Loading phase covariance data"
                try:
                    data=util.FITS.Read(fname)
                    covMatPhix=data[1]
                    self.Ax=data[3]
                    self.Bx=data[5]
                    self.AStartx=data[7]
                except:
                    print "Unable to load covariance data... generating"
                    traceback.print_exc()
                    covMatPhix=None
            if covMatPhix==None:
                print "Computation of the X phase covariance matrix"
                covMatPhix=self.computePhaseCovarianceMatrix(self.scrnXPxls+1,self.L0,self.pixScale,self.nbColToAdd,self.nbCol)
                print "Computation of the Ax and Bx matrixes"        
                self.Ax,self.Bx,self.AStartx=self.computeAandBmatrices(self.scrnXPxls+1,covMatPhix,self.nbColToAdd,self.nbCol)##we compute the A and B matrices
                util.FITS.Write(covMatPhix,fname)
                util.FITS.Write(self.Ax,fname,writeMode="a")
                util.FITS.Write(self.Bx,fname,writeMode="a")
                util.FITS.Write(self.AStartx,fname,writeMode="a")
            if self.saveCovMat:
                self.covMatPhix=covMatPhix
            else:
                del(covMatPhix)

            fname="scrn/infScrnData%d_%g_%g_%d_%d.fits"%(self.scrnYPxls+1,self.L0,self.pixScale,self.nbColToAdd,self.nbCol)
            covMatPhiy=None
            if os.path.exists(fname):
                print "Loading phase covariance data"
                try:
                    data=util.FITS.Read(fname)
                    covMatPhiy=data[1]
                    self.Ay=data[3]
                    self.By=data[5]
                    self.AStarty=data[7]
                except:
                    print "Unable to load covariance data... generating"
                    traceback.print_exc()
                    covMatPhiy=None
            if covMatPhiy==None:
                print "Computation of the Y phase covariance matrix"
                covMatPhiy=self.computePhaseCovarianceMatrix(self.scrnYPxls+1,self.L0,self.pixScale,self.nbColToAdd,self.nbCol)
                print "Computation of the Ay and By matrixes"        
                self.Ay,self.By,self.AStarty=self.computeAandBmatrices(self.scrnYPxls+1,covMatPhiy,self.nbColToAdd,self.nbCol)##we compute the A and B matrices
                util.FITS.Write(covMatPhiy,fname)
                util.FITS.Write(self.Ay,fname,writeMode="a")
                util.FITS.Write(self.By,fname,writeMode="a")
                util.FITS.Write(self.AStarty,fname,writeMode="a")
                
            if self.saveCovMat:
                self.covMatPhiy=covMatPhiy
            else:#save memory (this can be eg 1GB in size...)
                del(covMatPhiy)

            ##we compute the initial phase screen
            print "Computation of the initial phase screen"
##             if args.has_key("idstr"):
##                 idstr=args["idstr"]
##             else:
##                 idstr=None
            self.screen=computeInitialScreen(self.config,self.idstr[0])
            ##the field self.screen stores the na array
            ra.seed(self.seed)#reinitialise the seed for use here...
            
            ##this is the phase screen data.  Contains new cols then new rows
            #in dimension zero, and the data in dimension 1.
            if self.sendWholeScreen==0:
                self.outputData=na.zeros(self.dims,self.dataType)
                self.colOutput=self.outputData[:self.maxColAdd*(self.scrnYPxls+1)]#cmod.utils.arrayFromArray(self.outputData,(self.maxColAdd,self.scrnYPxls+1),self.dataType)
                self.colOutput.shape=(self.maxColAdd,self.scrnYPxls+1)
                self.rowOutput=self.outputData[self.maxColAdd*(self.scrnYPxls+1):]#cmod.utils.arrayFromArray(self.outputData[self.maxColAdd*(self.scrnYPxls+1):,],(self.maxRowAdd,self.scrnXPxls+1),self.dataType)
                self.rowOutput.shape=(self.maxRowAdd,self.scrnXPxls+1)
            else:
                self.outputData=self.screen
            print "TODO: infScrn - reduce memory requirements by only storing the part of the initial phasescreen that is required..."
            self.cmodInfo=None
            if self.useCmodule:
                nthreads=config.getVal("nthreads",default="all")
                if nthreads=="all":
                    nthreads=config.getVal("ncpu")
                self.randarr=numpy.zeros((self.scrnXPxls+1 if self.scrnXPxls>self.scrnYPxls else self.scrnYPxls+1,),numpy.float64)
                self.r0last=self.ro
                self.xsteplast=self.xstep
                self.ysteplast=self.ystep
                seed=0 if self.seed==None else self.seed
                self.cmodInfo=cmod.scrn.initialise(nthreads,self.ro,self.L0,self.scrnXPxls,self.scrnYPxls,self.sendWholeScreen,self.maxColAdd,self.maxRowAdd,self.colAdd,self.rowAdd,seed,self.screen,self.Ay,self.By,self.AStarty,self.Ax,self.Bx,self.AStartx,self.xstep,self.ystep,self.randarr,self.colOutput,self.rowOutput)

    def __del__(self):
        if self.cmodInfo!=None:
            cmod.scrn.free(self.cmodInfo)
        self.cmodInfo=None
        
    def computeR0(self):
        """Computes r0 as a function of time."""
        if self.r0Fun!=None:
            #r0Fun is a function that can return either a single float, or
            #an array of r0 (for turbulence moving at 0 or 90 degrees), of size 1+scrnXPxls or 1+scrnYPxls (depending if 90 or 0 degrees), or a tuple of 2 arrays, one for each direction for cases where not 90 degrees.  This allows phase dislocations etc.
            self.globR0=self.r0Fun(self.niter*self.tstep)
            if type(self.globR0)==type(()):
                self.ro=(self.globR0[0]*self.strLayerToPowMinusThreeOverFive,self.globR0[1]*self.strLayerToPowMinusThreeOverFive)
            else:
                self.ro=self.globR0*self.strLayerToPowMinusThreeOverFive
        if self.stepFun!=None:#add a step to the phase...
            if self.xstep==0 and self.ystep==0:#only if hasn't been added last time
                step=self.stepFun(self.niter*self.tstep)
                if type(step)==type(()):
                    self.ystep=step[1]
                    self.xstep=step[0]
                    print "got steps %s"%str(self.xstep)
                else:
                    self.xstep=self.ystep=step

    def computePhaseCovarianceMatrix(self,size,L0,pixScale,nbColToAdd=1,nbCol=2):
        """Computes the phase covariance matrix required to compute the <XXT>, <XZT> and <ZZT> matrices
        used to compute the A and B matrices
        cf equation 5 of Optics Express paper
        If nbCol==2 and nbColToAdd==1, the phase covariance will be a 3x3
        block matrix, with some blocks equal to others, ie
        (0,0 being bottom left)...
        a b c
        b c b
        c b a
        Each block equals its transpose, and so does the whole thing.
        Speed improvements by AGB.
        """
        ##we add one column at each iteration
        #nbColToAdd=1
        #self.nbColToAdd=nbColToAdd
        
        ##class characteristics
        #nbCol=self.nbCol#2
        dpix=size#agb changed from: self.dpix
        
        nc=nbCol+nbColToAdd
        ##we create a grid of points by ascending order in the vertical direction,
        ##from left to right
        nbPoints=nc*dpix
        
        ##numbers used for the phase structure matrix computation
        f0=1./L0;
        coeff=2*gamma(11./6)/(2**(5./6))/(na.pi**(8./3))
        coeff*=(24./5*gamma(6./5))**(5./6)
        
        ##phase variance for l0/r0=1 (scales as (l0/r0)**(5./3))
        sigma2=coeff/2.*gamma(5./6)/(2**(1./6))
        
        ##allocation of the covariance matrix
        covMatPhi=na.empty((nbPoints,nbPoints),na.float64,order="F")
        distMap=distanceMap((nc*2-1),dpix*2-1,-0.5,-0.5)*pixScale*2*na.pi*f0
        cPhi=distMap**(5./6)
        #same here
        cPhi*=kv(5./6,distMap)
        #util.FITS.Write(cPhi,"cPhi.fits")
        cPhi*=coeff/2.
        ##we create the covariance matrix
        for i in range(nbPoints):
            ##if not(i%dpix):
            ##print "Point",i,"sur",nbPoints-1
            ##Coordinates of the point we want to compute the phase covariance
            ic1=i%int(dpix)
            jc1=i/int(dpix)
            ##grid of distances centered on the point
            #distMap=distanceMap(dpix,nbCol+nbColToAdd,-dpix/2.+ic1,-nbCol/2.-nbColToAdd/2.+jc1)*self.pixScale
            ##we compute the phase covaraince function
            #cPhi=(2*na.pi*distMap*f0)**(5./6)
            #cPhi*=kv(5./6,2*na.pi*distMap*f0)
            #cPhi*=coeff/2.
            #cPhi[ic1,jc1]=sigma2 ##for the variance
            ##filling of the covariance matrix
            #print nc-1-jc1,nc*2-1-jc1,dpix-1-ic1,dpix*2-1-ic1,cPhi.shape,dpix,type(dpix)
            line=na.ravel(cPhi[nc-1-jc1:nc*2-1-jc1,dpix-1-ic1:dpix*2-1-ic1]).astype(covMatPhi.dtype)
            line[jc1*dpix+ic1]=sigma2#for the variance.
            covMatPhi[i,i:]=line[i:,]
            covMatPhi[i:,i]=line[i:,]
        return covMatPhi
    
    def computePhaseCovarianceMatrixOld(self,size):
        """Computes the phase covariance matrix required to compute the <XXT>, <XZT> and <ZZT> matrices
        used to compute the A and B matrices
        cf equation 5 of Optics Express paper
        If nbCol==2 and nbColToAdd==1, the phase covariance will be a 3x3
        block matrix, with some blocks equal to others, ie
        (0,0 being bottom left)...
        a b c
        b c b
        c b a
        Each block equals its transpose, and so does the whole thing.
        This is a slow one...(FA)
        """
        ##we add one column at each iteration
        nbColToAdd=1
        self.nbColToAdd=nbColToAdd
        
        ##class characteristics
        nbCol=self.nbCol
        dpix=size#agb changed from: self.dpix
        
        ##we create a grid of points by ascending order in the vertical direction,
        ##from left to right
        nbPoints=(nbCol+nbColToAdd)*dpix
        gridPoints=na.arange(nbPoints)
        gridPoints.shape=(nbCol+nbColToAdd,dpix)
        gridPoints=na.transpose(gridPoints)
        
        ##numbers used for the phase structure matrix computation
        f0=1./self.L0;
        coeff=2*gamma(11./6)/(2**(5./6))/(na.pi**(8./3))
        coeff*=(24./5*gamma(6./5))**(5./6)
        
        ##phase variance
        sigma2=coeff/2.*gamma(5./6)/(2**(1./6))
        
        ##allocation of the covariance matrix
        covMatPhi=na.empty((nbPoints,nbPoints),na.float64,order="F")
        
        ##we create the covariance matrix
        for i in range(nbPoints):
            ##if not(i%dpix):
            ##print "Point",i,"sur",nbPoints-1
            ##Coordinates of the point we want to compute the phase covariance
            ic1=i%dpix
            jc1=i/dpix
            ##grid of distances centered on the point
            distMap=distanceMap(dpix,nbCol+nbColToAdd,-dpix/2.+ic1,-nbCol/2.-nbColToAdd/2.+jc1)*self.pixScale
            ##we compute the phase covaraince function
            cPhi=(2*na.pi*distMap*f0)**(5./6)
            cPhi*=kv(5./6,2*na.pi*distMap*f0)
            cPhi*=coeff/2.
            cPhi[ic1,jc1]=sigma2 ##for the variance
            ##filling of the covariance matrix
            line=na.ravel(na.transpose(cPhi)).astype(covMatPhi.dtype)
            covMatPhi[i,i:]=line[i:,]
            covMatPhi[i:,i]=line[i:,]
        return covMatPhi
    
    def computeAandBmatrices(self,size,covMatPhi,nbColToAdd=1,nbCol=2):
        """Computes the A and B matrices from the phase covariance matrix
        """
        dpix=size#agb changed from: self.dpix
        N=dpix*(nbCol+nbColToAdd)
        M=dpix*nbColToAdd
        
        ##we declare the matrices with the fortran keyword to improve efficiency
        ##first matrix: ZZT
        #s=(self.nbCol*dpix,self.nbCol*dpix)#agb commented out on 070426
        #ZZT=na.empty(s,dtype=na.float64,order="F")
        ZZT=covMatPhi[:N-M,:N-M]
##         gist.window(0,wait=1)
##         gist.fma()
##         gist.pli(ZZT)
##         gist.fma()
##         raw_input()
        ##second matrix: XXT
        #s=(dpix,dpix)
        #XXT=na.empty(s,dtype=na.float64,order="F")
        XXT=covMatPhi[N-M:,N-M:]
        
        ##third matrix: XZT
        #s=(dpix,self.nbCol*dpix)
        #XZT=na.empty(s,dtype=na.float64,order="F")
        XZT=covMatPhi[N-M:,:N-M]
        
        ##fourth matrix: ZXT = transpose of XZT
        #s=(self.nbCol*dpix,dpix)
        #ZXT=na.empty(s,dtype=na.float64,order="F")
        ZXT=covMatPhi[:N-M,N-M:]

        ##we compute the inverse of ZZT
        t0=time.time()
        print "infScrn - doing cho_solve 0"
        #util.FITS.Write(ZZT,"zzt.fits")
        try:
            ZZT_inv=LA.cho_solve(LA.cho_factor(ZZT),na.identity(nbCol*dpix))
        except:
            print "cho_solve failed - trying inv... this sometimes happens if the matrix is too large... or is r0/pxl is too small."
            ZZT_inv=LA.inv(ZZT)
        
        ##we compute a Fortran contiguous matrix
        print "infScrn - doing matrix dot %g"%(time.time()-t0)
        matrixA=matrixdot(XZT,ZZT_inv)
        ##A=generalMatrixMatrixMultiply(XZT,ZZT_inv)
        #matrixA=na.empty(A.shape,dtype=na.float64,order="F")
        #matrixA[:,:]=A
        ##the A matrix computed there can be used to add new rows or new columns
        ##at the END of the phase screen
        ##I found that the matrix AStart used to add new rows or columns at the START is the vertical reciproc
        ##of the A matrix
        ##Ex : if A=(A2 A1 A0) then AStart=(A0 A1 A2)
        matrixAStart=na.empty(matrixA.shape,dtype=na.float64)#,order="F")
        ##we fill the matrix
        for col in range(nbCol):
            jstart=col*dpix
            jend=(col+1)*dpix
            ##print "jstart=%d jend=%d" % (jstart,jend)
            matrixAStart[:,jstart:jend]=matrixA[:,-(jstart+1):-(jend+1):-1][:,::-1]
        
        ##we compute now the B matrix
        print "infScrn - doing 2nd matrix dot %g"%(time.time()-t0)
        BBt=XXT-matrixdot(matrixA,ZXT)
        print "infScrn - doing mmx %g"%(time.time()-t0)
        ##BBt=XXT-generalMatrixMatrixMultiply(self.A,ZXT)
        print "infScrn - doing svd %g"%(time.time()-t0)
        u,w,vt=LA.svd(BBt)
        L=na.sqrt(w)
        matrixB=(u*L[na.newaxis,:]).copy()#multiplication of columns of U by diagonal elements
        #matrixB=na.empty(B.shape,dtype=na.float64,order="F")
        #matrixB[:,:]=B#this copy may be necessary, because C is c contiguous, not fortran.
        print "infScrn - done A and B matricees %g"%(time.time()-t0)
        return matrixA,matrixB,matrixAStart
    ##The next functions are the functions used to add new rows or columns at the beginning or the end of the phase screen
    def addNewColumnOnEnd(self,ro=None):
        """Updates the phase screen by adding a new column to the end of the phase screen and
        putting away the first one
        do the same than add addNewColumn
        """
        if ro is None:
            ro=self.ro
        if type(ro)==type(()):
            ro=ro[0]
        dpix=self.scrnYPxls+1#agb: changed from dpix
        ##we extract the last nbCol columns of the phase screen
        oldPhi=self.screen[:,-self.nbCol:]
        ##we put it into a single vector
        Z2=na.ravel(na.transpose(oldPhi)).astype(na.float64) ##Float64 for data precision
        ##multiplication by the A matrix (BLAS function)
        AZ=matrixdot(self.Ay,Z2)#agbc was Ay
        #txt="AZ: "
        #for i in range(AZ.shape[0]):
        #    txt+="%g, "%AZ[i]
        #print txt

        ##creation of random variables with the good variances
        coeffTurb=(self.L0/ro)**(5./6)
        #print "coeffTurb %s"%str(coeffTurb)
        rn=ra.randn(self.nbColToAdd*dpix)#By.shape==nbColToAdd*dpix...
        #txt="rn: "
        #for i in range(rn.shape[0]):
        #    txt+="%g, "%rn[i]
        #print txt
        rn=matrixdot(self.By,rn)*coeffTurb#agbc was By
        ##rn=generalMatrixVectorMultiply(self.B,rn)
        rn+=AZ ##vector storing the values of the last column
        if type(self.xstep)==na.ndarray or self.xstep!=0.:
            rn+=self.xstep
            self.xstep=0.
        ##we update the phase screen
        ##we first put out the first columns
        #print self.screen.shape,dpix,dpix-self.nbColToAdd
        self.screen[:,:-self.nbColToAdd]=self.screen[:,self.nbColToAdd:]#.copy()
        ##we then replace the last column
        #print self.screen.shape,rn.shape
        self.screen[:,-1]=rn.astype(self.screen.dtype)
        #txt=""
        #for i in range(self.screen.shape[1]):
        #    txt+="%g, "%self.screen[i,-1]
        #print txt
    def addNewColumnOnStart(self,ro=None):
        """Updates the phase screen by adding a new column to the start of the phase screen and
        putting away the last one
        """
        if ro is None:
            ro=self.ro
        if type(ro)==type(()):
            ro=ro[0]
        dpix=self.scrnYPxls+1#agb: changed from dpix
        ##we extract the first nbCol columns of the phase screen
        oldPhi=self.screen[:,:self.nbCol]
        ##we put it into a single vector
        Z2=na.ravel(na.transpose(oldPhi)).astype(na.float64) ##Float64 for data precision
        ##multiplication by the A matrix (BLAS function)
        AZ=matrixdot(self.AStarty,Z2)
        ##creation of random variables with the good variances
        coeffTurb=(self.L0/ro)**(5./6)
        rn=ra.randn(self.nbColToAdd*dpix)
        rn=matrixdot(self.By,rn)*coeffTurb
        ##rn=generalMatrixVectorMultiply(self.B,rn)
        rn+=AZ ##vector storing the values of the last column
        if type(self.xstep)==na.ndarray or self.xstep!=0.:
            rn+=self.xstep
            self.xstep=0.
        
        ##we update the phase screen
        ##we first put out the last columns
        self.screen[:,self.nbColToAdd:]=self.screen[:,:-self.nbColToAdd].copy()#the copy is needed, otherwise, overwrites itself giving incorrect results.
        ##we then replace the first column
        rn.shape=(self.screen.shape[0],self.nbColToAdd)
        self.screen[:,:self.nbColToAdd]=rn.astype(self.screen.dtype)

    def addNewRowOnEnd(self,ro=None):
        """Updates the phase screen by adding a new row to the end of the phase screen and
        putting away the first one
        """
        if ro is None:
            ro=self.ro
        if type(ro)==type(()):
            ro=ro[1]
        dpix=self.scrnXPxls+1#agb: changed from dpix
        
        ##we extract the last nbCol columns of the phase screen
        oldPhi=self.screen[-self.nbCol:,:]
        ##we put it into a single vector
        Z2=na.ravel(oldPhi).astype(na.float64) ##Float64 for data precision
        ##multiplication by the A matrix (BLAS function)
        AZ=matrixdot(self.Ax,Z2)
        ##AZ=generalMatrixVectorMultiply(self.A,Z2)
        ##creation of random variables with the good variances
        coeffTurb=(self.L0/ro)**(5./6)
        rn=ra.randn(self.nbColToAdd*dpix)
        rn=matrixdot(self.Bx,rn)*coeffTurb
        rn+=AZ ##vector storing the values of the last column
        if type(self.ystep)==na.ndarray or self.ystep!=0.:
            rn+=self.ystep
            self.ystep=0.

        ##we update the phase screen
        ##we first put out the first rows
        self.screen[:-self.nbColToAdd]=self.screen[self.nbColToAdd:,]#.copy()
        ##we then replace the last row
        self.screen[-1,:]=rn.astype(self.screen.dtype)

    def addNewRowOnStart(self,ro=None):
        """Updates the phase screen by adding a new row to the start of the phase screen and
        putting away the last one
        """
        if ro is None:
            ro=self.ro
        if type(ro)==type(()):
            ro=ro[1]
        dpix=self.scrnXPxls+1#agb: changed from dpix
        ##we extract the first nbCol columns of the phase screen
        oldPhi=self.screen[:self.nbCol,:]
        ##we put it into a single vector
        Z2=na.ravel(oldPhi).astype(na.float64) ##Float64 for data precision
        ##multiplication by the A matrix (BLAS function)
        AZ=matrixdot(self.AStartx,Z2)
        ##creation of random variables with the good variances
        coeffTurb=(self.L0/ro)**(5./6)
        rn=ra.randn(self.nbColToAdd*dpix)
        rn=matrixdot(self.Bx,rn)*coeffTurb
        ##rn=generalMatrixVectorMultiply(self.B,rn)
        rn+=AZ ##vector storing the values of the last column
        if type(self.ystep)==na.ndarray or self.ystep!=0.:
            rn+=self.ystep
            self.ystep=0.
            
        ##we update the phase screen
        ##we first put out the last rows
        self.screen[self.nbColToAdd:,]=self.screen[:-self.nbColToAdd,]#.copy()
        ##we then replace the first column
        self.screen[:self.nbColToAdd,:]=rn.astype(self.screen.dtype)

    def addNewData(self,ro=None):
        """Updates the phase screen by adding new rows or columns
        as specified by the wind direction
        """
        self.computeR0()
        if ro is None:
            ro=self.ro
        if self.useCmodule:
            if ro!=self.r0last:
                if type(ro)==type(()):
                    cmod.scrn.update(self.cmodInfo,2,ro[0])
                    cmod.scrn.update(self.cmodInfo,3,ro[1])
                else:
                    cmod.scrn.update(self.cmodInfo,1,ro)
                self.r0last=ro
            if self.xstep!=self.xsteplast:
                cmod.scrn.update(self.cmodInfo,4,self.xstep)
                self.xsteplast=self.xstep
            if self.ystep!=self.ysteplast:
                cmod.scrn.update(self.cmodInfo,5,self.ystep)
                self.ysteplast=self.ystep
            #print self.cmodInfo,self.cmodInfo.flags,self.cmodInfo.size,self.cmodInfo.itemsize,self.cmodInfo.dtype.char
            #randn=ra.randn(self.scrnXPxls+1 if self.scrnXPxls>self.scrnYPxls else self.scrnYPxls+1)
            #cmod.scrn.update(self.cmodInfo,6,randn)
            cmod.scrn.run(self.cmodInfo)
        else:
            nrem,nadd,interppos=self.newCols.next()
            nadd=int(nadd)
            # if colAdd<zero, we are adding new columns on the right of the array.
            # If rowAdd<zero, we are adding new rows at the bottom of the array (note, that if plotting in Gist, this is the top of the array).
            if self.colAdd<0:
                for i in range(nadd):
                    self.addNewColumnOnEnd(ro)
            else:
                for i in range(nadd):
                    self.addNewColumnOnStart(ro)
            nrem,nadd,interppos=self.newRows.next()
            nadd=int(nadd)
            if self.rowAdd<0:
                for i in range(nadd):
                    self.addNewRowOnEnd(ro)
            else:
                for i in range(nadd):
                    self.addNewRowOnStart(ro)
    def prepareOutput(self):
        """If not sending whole screen, copies the parts to be sent..."""
        if self.sendWholeScreen==0:
            if self.colAdd<0:
                self.colOutput[:,:]=na.transpose(self.screen[:,-self.maxColAdd:])
            else:
                self.colOutput[:,:]=na.transpose(self.screen[:,:self.maxColAdd])
            if self.rowAdd<0:
                self.rowOutput[:,:]=self.screen[-self.maxRowAdd:,]
            else:
                self.rowOutput[:,:]=self.screen[:self.maxRowAdd]
    
    def generateNext(self,msg=None): #FA screen main loop
        """
        This function is called when it is okay to produce the next iteration
        of the simulation.
        Not expecting any msgs
        """
        t1=time.time()
        ##we first do a loop over the pixel offset
        if self.generate==1: ##Other modules call this module
            self.dataValid=1 ##Data is OK
            self.niter+=1
            self.addNewData()
            self.prepareOutput()
        else:
            self.dataValid=0 ##No new data is required
        self.generateNextTime=time.time()-t1
        

    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr==None:
            id=""
        else:
            id=" (%s)"%self.idstr
        return """<plot title="Phase screen%s" cmd="data=%s.screen" ret="data" type="pylab" when="rpt" palette="gray"/>"""%(id,objname)
	
    def getParams(self):
        """parameters required for this module, in the form of {"paramName":defaultValue,...}
        These params can then be placed in the config file... if not set by the
        user, the param should still be in config file as default value for
        future reference purposes.
        """
        #This is a working example.  Please feel free to change the parameters
        #required. (if you do, also change the config.getParam() calls too).
        paramList=[]
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="dataType",typ="eval",val="this.globals.fpDataType",comment="Array numpy data type"))
        paramList.append(base.dataType.dataType(description="l0",typ="f",val="30.0",comment="TODO: turbulence outer scale"))
        paramList.append(base.dataType.dataType(description="r0",typ="f",val="0.15",comment="TODO: Fried parameter (m)"))
        #paramList.append(base.dataType.dataType(description="strength",typ="eval",val="this.infScrn.strenghts['0m']",comment="TODO: layer strength (a float, initially a dictionary holding all strengths)"))
        #paramList.append(base.dataType.dataType(description="windDirection",typ="eval",val="this.infAtmos.windDirection['0m']",comment="TODO: Wind direction (degrees, going from -180 to 180) - include as dict in infAtmos module, with keys equal to layer heights, and then select the correct one for each infScrn module."))
        #paramList.append(base.dataType.dataType(description="vWind",typ="eval",val="this.infAtmos.vWind['0m']",comment="TODO: Wind velocity (m/s) - include as dict in infAtmos module, with keys equal to layer heights, and then select the correct one for each infScrn module."))
        #paramList.append(base.dataType.dataType(description="offset",typ="i",val="1",comment="wind offset in pixels (integer)"))
        paramList.append(base.dataType.dataType(description="nbCol",typ="i",val="2",comment="Number of columns or rows used to create a new one"))
        #paramList.append(base.dataType.dataType(description="initPhsSeed",typ="eval",val="None",comment="Random number generator seed"))
        #paramList.append(base.dataType.dataType(description="scrnSize",typ="code",val="from science.infScrn import computeScrnSize;scrnSize=computeScrnSize(this.infAtmos.sourceThetaDict,this.infAtmos.sourcePhiDict,this.globals.ntel,this.globals.npup,this.globals.telDiam,this.infAtmos.altitude)",comment="Screen sizes"))
        #paramList.append(base.dataType.dataType(description="scrnXPxls",typ="eval",val="this.infScrn.scrnSize['0m'][0]",comment="TODO: select the screen size for this layer."))        
        #paramList.append(base.dataType.dataType(description="scrnYPxls",typ="eval",val="this.infScrn.scrnSize['0m'][1]",comment="TODO: select the screen size for this layer."))        
        paramList.append(base.dataType.dataType(description="tstep",typ="f",val="0.005",comment="TODO: timestep."))        
        #paramList.append(base.dataType.dataType(description="degRad",typ="eval",val="2*Numeric.pi/360.",comment="degrees to radians."))        
        paramList.append(base.dataType.dataType(description="stepFunction",typ="eval",val="None",comment="For modelling phasesteps."))        
        paramList.append(base.dataType.dataType(description="r0Function",typ="eval",val="None",comment="R0 as a function of time."))        
        paramList.append(base.dataType.dataType(description="saveInfPhaseCovMatrix",typ="i",val="0",comment="Save the inf phase covariance matrix."))        
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))

        return paramList
	
    def getInputType(self):
        """Returns the input needed by this module.
        Should return an instance of dataType, or a list of dataType objects,
        or None if no input is required."""
        return None
	
    def getOutputType(self):
        """Returns the output given by this module.
        Should return an instance of dataType, or a list of dataType objects,
        or None if no output is given."""
        return base.dataType.dataType("phasescreen",na.ndarray,"*","*")
	
    def getInitArgs(self):
        """return a dictionary of parameters that can be passed to init.
        Dictionary key is the parameter name, and dictionary value is a tuple
        of (description, default value).  If there is no default value, then
        the dictionary value is a tuple of (description,).
        """
        return {"config":("Config object",),"forGUISetup":("initialisation type",0)}

if __name__=="__main__":
    import base.readConfig
    config=base.readConfig.AOXml("params.xml")
    scrn1=Mkscrns(None,config,args={"idstr":1}) ##Create the object
    mkscrns.next()
    mkscrns.next()
