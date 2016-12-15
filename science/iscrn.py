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
#todo: could be made more efficient by adding all new columns/rows at once,
#rather than shifting, adding, shifting, adding... this would save on shifts.
#Not sure how much of a speedup this would give.
#todo: think about whether its worth passing the whole phasescreen, rather than
#just the new rows each time - if the infAtmos object is on the same node, this
#will save computation in infAtmos, and take no time to pass.  However, if on
#a different node, may take longer.  Maybe, we need some option for this.
#Maybe outputData should contain both, which can then be selected.  Or maybe
#infAtmos should look at parent.screen to see if it exists, and if so, use it.
#Yes, thats probably the most elegant way...

##Infinite phase screen module, with wrap implemented.
"""iscrn.py:  Infinite phase screen generation using
Francois Assemat's infinite phase screens method
for phase screen computation

cf Assemat et al., Optics Express, Vol. 14, Issue 3, pp. 988-999 (February 2006)
Things to think about:

Creates new rows of atmospheric phase each iteration.  No
interpolation is done.  The number of rows added depends on
wind velocity.  If a non-integer number is required each
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
import cmod.iscrn

## for fast MVM:
import util.dot as quick # will replace numpy.dot when I implement it for double

matrixdot=quick.dot#matrix.dot

    
# def computeScrnSize(thetas,phis,ntel,npup,telDiam,altitude,zenith):
#     """Zenith in degrees"""
    
#     zenith=0.
#     arcsecrad=2*na.pi/360./3600.
#     degrad=2*na.pi/360.
#     scrnSize={}
#     for altkey in altitude.keys():
#         xposlist=[]
#         yposlist=[]
#         layerAlt=altitude[altkey]
#         for key in thetas.keys():
#             xposlist.append(layerAlt/numpy.cos(zenith*degrad)*na.fabs(na.tan(thetas[key]*arcsecrad)*na.cos(phis[key]*degrad)))
#             yposlist.append(layerAlt*na.fabs(na.tan(thetas[key]*arcsecrad)*na.sin(phis[key]*degrad)))
#         maxx=max(xposlist)
#         minx=min(xposlist)
#         maxy=max(yposlist)
#         miny=min(yposlist)
#         scrnXPxls=int(na.ceil(maxx*ntel/telDiam+npup+na.ceil(minx*ntel/telDiam))+1)
#         scrnYPxls=int(na.ceil(maxy*ntel/telDiam+npup+na.ceil(miny*ntel/telDiam))+1)
#         scrnSize[altkey]=(scrnXPxls,scrnYPxls)
#     print "scrnSize: %s"%str(scrnSize)
#     return scrnSize

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
        searchOrder=["iscrn_"+idstr,"iscrn","globals"]
    else:
        searchOrder=["iscrn","globals"]
    config.setSearchOrder(searchOrder)
    Dtel=config.getVal("telDiam")#diameter in m.
    dpix=config.getVal("npup")#diameter in pixels.
    atmosGeom=config.getVal("atmosGeom",default=None,raiseerror=1)
    windDirection=atmosGeom.layerWind(idstr)
    scrnXPxls=atmosGeom.getScrnXPxls(idstr,rotateDirections=1)
    scrnYPxls=atmosGeom.getScrnYPxls(idstr,rotateDirections=1)
    seed=atmosGeom.layerInitSeed(idstr)
    strLayer=atmosGeom.layerStrength(idstr)
    vWind=atmosGeom.layerSpeed(idstr)
    globR0=atmosGeom.r0
    L0=atmosGeom.l0
    print "Generating phasescreen %s size %s"%(str(idstr),str((scrnXPxls,scrnYPxls)))
    tstep=config.getVal("tstep")
    scrnDir=config.getVal("scrnDir",default="scrn")
    natype="d"#config.getVal("dataType")
    phaseArray=makeInitialScreen(dpix,Dtel,L0,scrnXPxls,scrnYPxls,seed,tstep,globR0,strLayer,natype,windDirection,vWind,scrnDir=scrnDir)
    ##we return the array corresponding to the pupil
    config.setSearchOrder(so)
    return phaseArray

def makeScrnQuick(npup,telDiam,l0=30.,r0=0.2,seed=0,scrnDir="scrn"):
    """For users on the command line, to generate a quick phase screen"""
    raise Exception("Need to sort this - clipping/oversizing no longer necessary")
    scrn=makeInitialScreen(dpix=npup+1,Dtel=telDiam,L0=l0,globR0=r0,seed=seed,scrnDir=scrnDir)[1:-1,1:-1]
    return scrn
def makeInitialScreen(dpix=1024,Dtel=42.,L0=30.,scrnXPxls=None,scrnYPxls=None,seed=0,tstep=0.05,globR0=0.2,strLayer=1.,natype=na.float64,windDirection=0.,vWind=10.,scrnDir="scrn"):
    """dpix is the telescope aperture diameter, and dtel is tel diameter.
    The actual number of pixels used is scrnXPxls x scrnYPxls.
    """

    if scrnXPxls==None:
        scrnXPxls=dpix
    if scrnYPxls==None:
        scrnYPxls=dpix
    if seed!=None and scrnDir!=None:
        try:
            fname=os.path.join(scrnDir,"iscrn%dD%gL%gx%dy%ds%dt%gr%gs%gd%gv%g.fits"%(dpix,Dtel,L0,scrnXPxls,scrnYPxls,seed,tstep,globR0,strLayer,windDirection,vWind))#windDirection needed in name because it determines scrnXPxls (in atmosGeom computation).
        except:
            fname=None
            print "Could not make filename for screen - so generating screen..."
    else:#random seed
        fname=None
    if fname!=None and os.path.exists(fname):
        print "Loading existing screen: %s"%fname
        try:
            phaseArray=util.FITS.Read(fname)[1]
            return phaseArray
        except:
            print "Unable to load file %s - generating"%fname
            fname=None

    scrnPxls=max(scrnXPxls,scrnYPxls)
    pixScale=Dtel/float(dpix)
    r0=globR0*(strLayer**(-3./5.))##we compute the r0 in the considered layer
    #colAdd=-vWind*na.cos(windDirection*na.pi/180)/pixScale*tstep#number of pixels to step each iteration (as float).
    #rowAdd=-vWind*na.sin(windDirection*na.pi/180)/pixScale*tstep#number of pixels to step each iteration (as float).
    #colAdd=0.#we don't add columns any more
    #rowAdd=-vWind/pixScale*tstep

    if seed==None:
        print "ERROR (possibly): computeInitialScreen - seed is None, so timer will be used, meaning that the initial screen cannot be replicated, so if both iscrn and infAtmos try to create, you will get a bad phasescreen.  If you wish to use a random seed, use int(time.time()) in the parameter file - though this will only work if all running in the same process."
    if L0>=Dtel:
        scrnSize=2*L0
    elif Dtel>=2*L0:
        scrnSize=Dtel
    else:
        scrnSize=2*L0
    ##size in pixels of the required phase screen
    nfft=int(na.round(dpix*1.*scrnSize/Dtel))#agb: scrnPxls was dpix
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
    f*=na.sqrt(fact*((scrnSize/r0)**(5./3.))*(nfft*1.)**2.) ## Von Karman phase power spectrum : 0.023/ro*(f^2+fo^2)^(-11/6)

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

    #phaseArray=na.zeros((scrnYPxls+1,scrnXPxls+1),natype)#was dpix agbc swapped
    #phaseArray[1:,:-1]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    phaseArray=phi[:scrnYPxls,:scrnXPxls].copy()
    # if colAdd<0:
    #     if rowAdd<0:
    #         phaseArray[1:,1:]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    #     else:
    #         phaseArray[:-1,1:]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    # else:
    #     if rowAdd<0:
    #         phaseArray[1:,:-1]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    #     else:
    #         phaseArray[:-1,:-1]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    if fname!=None and scrnDir!=None:
        if not os.path.exists(scrnDir):
            os.mkdir(scrnDir)
        util.FITS.Write(phaseArray,fname)
    return phaseArray

class This:
    pass

class iscrn(base.aobase.aobase):
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
        ##Physical size in meters per pixel
        self.pixScale=self.Dtel/float(self.dpix)
        ##Precision of the output screen array
        self.dataType="d"#config.getVal("dataType")
        if self.dataType not in ["d",na.float64]:
            print "WARNING - iscrn - it is not known whether statistics for iscrn will be right if not Float64.  Would be good if someone could test this... ie is it true Kolmogorov/von Karman"
        self.tstep=config.getVal("tstep")
        self.useCmodule=config.getVal("useCmodule",default=1)
        self.sendWholeScreen=self.config.getVal("sendWholeScreen",default=0)
        self.atmosGeom=config.getVal("atmosGeom")
        self.L0=self.atmosGeom.l0
        self.globR0=self.atmosGeom.r0
        #Is this object doing 1 layer, or several/many layers?
        self.r0Fun=config.getVal("r0Function",default=None,raiseerror=0)#function giving r0 as function of time.
        self.layerList=self.config.getVal("layerList",{}).get(self.idstr[0],None)
        if self.layerList==None:
            self.layerList=[self.idstr[0]]
        self.thisObjDict={}
        dim=0
        for id in self.layerList:#for each layer that I'm generating here...
            this=This()
            this.insertPos=0
            self.thisObjDict[id]=this
            # Total size of phasescreen (will be larger than dpix if there are off axis sources, and will also depend on height of the layer).
            this.scrnXPxls=self.atmosGeom.getScrnXPxls(id,rotateDirections=1)
            this.scrnYPxls=self.atmosGeom.getScrnYPxls(id,rotateDirections=1)
            this.windDirection=self.atmosGeom.layerWind(id)
            this.vWind=self.atmosGeom.layerSpeed(id)
            this.strLayer=self.atmosGeom.layerStrength(id)
            this.strLayerToPowMinusThreeOverFive=this.strLayer**(-3./5)
            this.windDirRad=this.windDirection*self.degRad
            # if colAdd<zero, we are adding now columns on the right of the array.
            # we are adding new rows at the bottom of the array (note, that if plotting in Gist, this is the top of the array).
            this.rowAdd=this.vWind/self.pixScale*self.tstep#number of pixels to step each iteration (as float).
            #This then tells us the size of the output array -
            this.maxRowAdd=int(na.ceil(na.fabs(this.rowAdd)))
            if self.sendWholeScreen==0:
                dim+=(this.maxRowAdd*(this.scrnXPxls))#was dpix
            else:
                dim+=(this.scrnYPxls*this.scrnXPxls)
            this.stepFun=config.getVal("stepFunction",default=None,raiseerror=0)#function giving phase steps... as function of time.
            if type(this.stepFun)==type({}):
                this.stepFun=this.stepFun[id]
            this.ystep=0.
            # Create the objects used to tell us how many extra cols/rows to add
            this.newRows=util.getNewCols.getNewCols(na.fabs(this.rowAdd))

        ##we compute the ro in the considered layer
        self.computeR0()

        self.dims=(dim,)
        ##Test for the future config GUI
        if forGUISetup==1:
            self.outputData=(self.dims,self.dataType)##specify output for GUI
        else: ##we continue the initialisation
            ##we extract the outer scale, the global Fried Parameter
            ##and the turbulence strength in the layer
            
            ##Number of previous rows or columns used to compute the new one (default value : 2)
            self.nbCol=2#config.getVal("nbCol",default=2)
            ##we extract the random number generator seed (for the additional rows/cols. At the moment, same for each layer used here.
            self.seed=config.getVal("seed",default=None,raiseerror=0)
            self.keepCovMat=config.getVal("keepInfPhaseCovMatrix",default=0)
            scrnDir=config.getVal("scrnDir",default="scrn")
            if not os.path.exists(scrnDir):
                os.makedirs(scrnDir)
            ##we go now through the creation of the required matrices
            ##we compute first the phase covariance matrices
            for id in self.layerList:
                this=self.thisObjDict[id]
                fname=os.path.join(scrnDir,"iscrnData%d_%g_%g.fits"%(this.scrnXPxls,self.L0,self.pixScale))
                covMatPhix=None
                if os.path.exists(fname):
                    print "Loading phase covariance data"
                    try:
                        data=util.FITS.Read(fname)
                        covMatPhix=data[1]
                        this.Ax=data[3]
                        this.Bx=data[5]
                        this.AStartx=data[7]
                    except:
                        print "Unable to load covariance data... generating"
                        traceback.print_exc()
                        covMatPhix=None
                if covMatPhix==None:
                    print "Computation of the X phase covariance matrix"
                    covMatPhix=self.computePhaseCovarianceMatrix(this.scrnXPxls,self.L0,self.pixScale,self.nbColToAdd,self.nbCol)
                    print "Computation of the Ax and Bx matrixes"        
                    this.Ax,this.Bx,this.AStartx=self.computeAandBmatrices(this.scrnXPxls,covMatPhix,self.nbColToAdd,self.nbCol)##we compute the A and B matrices
                    try:
                        util.FITS.Write(covMatPhix,fname)
                        util.FITS.Write(this.Ax,fname,writeMode="a")
                        util.FITS.Write(this.Bx,fname,writeMode="a")
                        util.FITS.Write(this.AStartx,fname,writeMode="a")
                    except:
                        print "Failed to write covariance matrices - disk full?  Continuing anyway..."
                if self.keepCovMat:
                    this.covMatPhix=covMatPhix
                else:
                    del(covMatPhix)
                print "Computation of the initial phase screen"
                this.screen=computeInitialScreen(self.config,id)


            ##we compute the initial phase screen

            ##the field self.screen stores the na array
            ra.seed(self.seed)#reinitialise the seed for use here...
            
            ##this is the phase screen data.  Contains new cols then new rows
            #in dimension zero, and the data in dimension 1.
            if self.sendWholeScreen==0:
                self.outputData=na.zeros(self.dims,self.dataType)
                #self.colOutput=self.outputData[:self.maxColAdd*(self.scrnYPxls+1)]#cmod.utils.arrayFromArray(self.outputData,(self.maxColAdd,self.scrnYPxls+1),self.dataType)
                #self.colOutput.shape=(self.maxColAdd,self.scrnYPxls+1)
                pos=0
                for id in self.layerList:
                    this=self.thisObjDict[id]
                    this.rowOutput=self.outputData[pos:pos+this.scrnXPxls*this.maxRowAdd]
                    this.rowOutput.shape=this.maxRowAdd,this.scrnXPxls
                    pos+=this.scrnXPxls*this.maxRowAdd
            else:
                self.outputData=numpy.zeros(self.dims,self.dataType)#self.screen
                pos=0
                for id in self.layerList:
                    this=self.thisObjDict[id]
                    tmp=self.outputData[pos:pos+this.scrnXPxls*this.scrnYPxls]
                    tmp.shape=this.scrnYPxls,this.scrnXPxls
                    tmp[:]=this.screen#copy the phase screen into the output array.
                    this.screen=tmp
                    pos+=this.scrnXPxls*this.scrnYPxls
            if self.useCmodule:
                nthreads=config.getVal("nthreads",default="all")
                if nthreads=="all":
                    nthreads=config.getVal("ncpu")
                for id in self.layerList:
                    this=self.thisObjDict[id]
                    this.randarr=numpy.zeros((this.scrnXPxls,),numpy.float64)# if self.scrnXPxls>self.scrnYPxls else self.scrnYPxls+1,),numpy.float64)
                    this.r0last=this.r0
                    this.ysteplast=this.ystep
                    seed=0 if self.seed==None else self.seed
                    this.cmodInfo=cmod.iscrn.initialise(nthreads,this.r0,self.L0,this.scrnXPxls,this.scrnYPxls,self.sendWholeScreen,this.maxRowAdd,this.rowAdd,seed,this.screen,this.Ax,this.Bx,this.AStartx,this.ystep,this.randarr,this.rowOutput)

    def __del__(self):
        for id in self.layerList:
            this=self.thisObjDict[id]
            if hasattr(this,"cmodInfo") and this.cmodInfo is not None:
                cmod.iscrn.free(this.cmodInfo)
            this.cmodInfo=None
        
    def computeR0(self):
        """Computes r0 as a function of time."""
        if self.r0Fun!=None:
            #r0Fun is a function that can return either a single float, or
            #an array of r0 (for turbulence moving up), of size scrnXPxls 
            self.globR0=self.r0Fun(self.niter*self.tstep)

        for id in self.layerList:
            this=self.thisObjDict[id]
            this.r0=self.globR0*this.strLayerToPowMinusThreeOverFive
            if this.stepFun!=None:#add a step to the phase...
                if this.ystep==0:#only if hasn't been added last time
                    this.ystep=this.stepFun(self.niter*self.tstep)

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
            ##Coordinates of the point we want to compute the phase covariance
            ic1=i%int(dpix)
            jc1=i//int(dpix)
            line=na.ravel(cPhi[nc-1-jc1:nc*2-1-jc1,dpix-1-ic1:dpix*2-1-ic1]).astype(covMatPhi.dtype)
            line[jc1*dpix+ic1]=sigma2#for the variance.
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
        ZZT=covMatPhi[:N-M,:N-M]
        XXT=covMatPhi[N-M:,N-M:]
        
        ##third matrix: XZT
        XZT=covMatPhi[N-M:,:N-M]
        
        ##fourth matrix: ZXT = transpose of XZT
        ZXT=covMatPhi[:N-M,N-M:]

        ##we compute the inverse of ZZT
        t0=time.time()
        print "iscrn - doing cho_solve 0"
        try:
            ZZT_inv=LA.cho_solve(LA.cho_factor(ZZT),na.identity(nbCol*dpix))
        except:
            print "cho_solve failed - trying inv... this sometimes happens if the matrix is too large... or is r0/pxl is too small."
            ZZT_inv=LA.inv(ZZT)
        
        ##we compute a Fortran contiguous matrix
        print "iscrn - doing matrix dot %g"%(time.time()-t0)
        matrixA=matrixdot(XZT,ZZT_inv)
        matrixAStart=na.empty(matrixA.shape,dtype=na.float64)#,order="F")
        ##we fill the matrix
        for col in range(nbCol):
            jstart=col*dpix
            jend=(col+1)*dpix
            ##print "jstart=%d jend=%d" % (jstart,jend)
            matrixAStart[:,jstart:jend]=matrixA[:,-(jstart+1):-(jend+1):-1][:,::-1]
        
        ##we compute now the B matrix
        print "iscrn - doing 2nd matrix dot %g"%(time.time()-t0)
        BBt=XXT-matrixdot(matrixA,ZXT)
        print "iscrn - doing svd %g"%(time.time()-t0)
        u,w,vt=LA.svd(BBt)
        L=na.sqrt(w)
        matrixB=(u*L[na.newaxis,:]).copy()#multiplication of columns of U by diagonal elements
        print "iscrn - done A and B matricees %g"%(time.time()-t0)
        return matrixA,matrixB,matrixAStart

    ##The next functions are the functions used to add new rows at the end of the phase screen

    def addNewRow(self,this,r0):
        """Updates the phase screen by adding a new row to the end (actually, at the current update position) of the phase screen, overwriting whatever is there.        """
        if r0 is None:
            r0=this.r0
        dpix=this.scrnXPxls#agb: changed from dpix
        ip=this.insertPos#position at which to insert new phase.
        ##we extract the last nbCol columns of the phase screen, dotting them with Ax to give AZ.  
        AZ=None
        for i in range(self.nbCol):#usually 2.
            indx=ip-self.nbCol+i
            if indx<0:#wrap around.
                indx+=this.scrnYPxls
            oldPhi=this.screen[indx]
            if AZ==None:
                AZ=matrixdot(this.Ax[:,i*this.scrnXPxls:(i+1)*this.scrnXPxls],oldPhi)
            else:
                AZ+=matrixdot(this.Ax[:,i*this.scrnXPxls:(i+1)*this.scrnXPxls],oldPhi)

        #oldPhi=self.screen[-self.nbCol:,:]
        ##we put it into a single vector
        #Z2=na.ravel(oldPhi).astype(na.float64) ##Float64 for data precision
        ##multiplication by the A matrix (BLAS function)
        #AZ=matrixdot(self.Ax,Z2)
        ##creation of random variables with the good variances
        coeffTurb=(self.L0/r0)**(5./6)
        rn=ra.randn(self.nbColToAdd*dpix)
        rn=matrixdot(this.Bx,rn)*coeffTurb
        rn+=AZ ##vector storing the values of the last column
        if type(this.ystep)==na.ndarray or this.ystep!=0.:
            rn+=this.ystep
            this.ystep=0.

        ##we update the phase screen
        ##we first put out the first rows
        #self.screen[:-self.nbColToAdd]=self.screen[self.nbColToAdd:,]#.copy()
        ##we then replace the active row
        this.screen[this.insertPos]=rn
        this.insertPos+=1
        if this.insertPos>=this.scrnYPxls:#wrap.
            this.insertPos=0


    def addNewData(self,userr0=None):
        """Updates the phase screen by adding new rows or columns
        as specified by the wind direction
        """
        self.computeR0()
        if self.useCmodule:
            for id in self.layerList:
                this=self.thisObjDict[id]
                if userr0 is None:
                    r0=this.r0
                else:
                    r0=userr0
                if r0!=this.r0last:
                    cmod.iscrn.update(this.cmodInfo,1,r0)
                    this.r0last=r0
                if this.ystep!=this.ysteplast:
                    cmod.iscrn.update(this.cmodInfo,5,this.ystep)
                    this.ysteplast=this.ystep
                this.insertPos=cmod.iscrn.run(this.cmodInfo)
        else:
            for id in self.layerList:
                this=self.thisObjDict[id]
                if userr0 is None:
                    r0=this.r0
                else:
                    r0=userr0
                nrem,nadd,interppos=this.newRows.next()
                nadd=int(nadd)
                for i in range(nadd):
                    self.addNewRow(this,r0)
            self.prepareOutput()

    def prepareOutput(self):
        """If not sending whole screen, copies the parts to be sent..."""
        if self.sendWholeScreen==0:
            for id in self.layerList:
                this=self.thisObjDict[id]
                for i in range(self.maxRowAdd):
                    pos=this.insertPos-self.maxRowAdd+i
                    if pos<0:
                        pos+=this.scrnYPxls#wrap.
                    this.rowOutput[i]=this.screen[pos]

    
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
        else:
            self.dataValid=0 ##No new data is required
        self.generateNextTime=time.time()-t1
    def unwrapPhase(self,id):
        this=self.thisObjDict[id]
        return numpy.concatenate([this.screen[this.insertPos:],this.screen[:this.insertPos]])

    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        txt=""
        for id in self.layerList:
            txt+="""<plot title="Raw Phase screen%s" cmd="data=%s.thisObjDict['%s'].screen" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,id)
            txt+="""<plot title="Unwrapped phase screen%s" cmd="data=%s.unwrapPhase('%s')" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,id)
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
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="dataType",typ="eval",val="this.globals.fpDataType",comment="Array numpy data type"))
        paramList.append(base.dataType.dataType(description="l0",typ="f",val="30.0",comment="TODO: turbulence outer scale"))
        paramList.append(base.dataType.dataType(description="r0",typ="f",val="0.15",comment="TODO: Fried parameter (m)"))
        #paramList.append(base.dataType.dataType(description="strength",typ="eval",val="this.iscrn.strenghts['0m']",comment="TODO: layer strength (a float, initially a dictionary holding all strengths)"))
        #paramList.append(base.dataType.dataType(description="windDirection",typ="eval",val="this.infAtmos.windDirection['0m']",comment="TODO: Wind direction (degrees, going from -180 to 180) - include as dict in infAtmos module, with keys equal to layer heights, and then select the correct one for each iscrn module."))
        #paramList.append(base.dataType.dataType(description="vWind",typ="eval",val="this.infAtmos.vWind['0m']",comment="TODO: Wind velocity (m/s) - include as dict in infAtmos module, with keys equal to layer heights, and then select the correct one for each iscrn module."))
        #paramList.append(base.dataType.dataType(description="offset",typ="i",val="1",comment="wind offset in pixels (integer)"))
        paramList.append(base.dataType.dataType(description="nbCol",typ="i",val="2",comment="Number of columns or rows used to create a new one"))
        #paramList.append(base.dataType.dataType(description="initPhsSeed",typ="eval",val="None",comment="Random number generator seed"))
        #paramList.append(base.dataType.dataType(description="scrnSize",typ="code",val="from science.iscrn import computeScrnSize;scrnSize=computeScrnSize(this.infAtmos.sourceThetaDict,this.infAtmos.sourcePhiDict,this.globals.ntel,this.globals.npup,this.globals.telDiam,this.infAtmos.altitude)",comment="Screen sizes"))
        #paramList.append(base.dataType.dataType(description="scrnXPxls",typ="eval",val="this.iscrn.scrnSize['0m'][0]",comment="TODO: select the screen size for this layer."))        
        #paramList.append(base.dataType.dataType(description="scrnYPxls",typ="eval",val="this.iscrn.scrnSize['0m'][1]",comment="TODO: select the screen size for this layer."))        
        paramList.append(base.dataType.dataType(description="tstep",typ="f",val="0.005",comment="TODO: timestep."))        
        #paramList.append(base.dataType.dataType(description="degRad",typ="eval",val="2*Numeric.pi/360.",comment="degrees to radians."))        
        paramList.append(base.dataType.dataType(description="stepFunction",typ="eval",val="None",comment="For modelling phasesteps."))        
        paramList.append(base.dataType.dataType(description="r0Function",typ="eval",val="None",comment="R0 as a function of time."))        
        paramList.append(base.dataType.dataType(description="saveInfPhaseCovMatrix",typ="i",val="0",comment="Save the inf phase covariance matrix."))        
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by iscrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))

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
