import numpy
import cmod.imgnoise
import util.atmos 
import util.tel
import util.FITS
import base.readConfig
import scipy.interpolate
#import sor30
#import plwf

def make_circle(r):
    t = numpy.arange(0, numpy.pi * 2.0, 0.01)
    t = t.reshape((len(t), 1))
    x = r * numpy.cos(t)
    y = r * numpy.sin(t)
    return numpy.hstack((x, y))
def makeDoughnut(outer,inner,xoff,yoff):
    import matplotlib.path as mpath
    vert=make_circle(outer)
    codes=numpy.ones(len(vert),dtype=mpath.Path.code_type)*mpath.Path.LINETO
    codes[0]=mpath.Path.MOVETO
    if inner!=None and inner!=0:
        inn=make_circle(inner)
        vert=numpy.concatenate((vert,inn[::-1]))
        codes=numpy.concatenate((codes,codes))
    vert[:,0]+=xoff
    vert[:,1]+=yoff
    path=mpath.Path(vert,codes)
    return path

def displayGSOverlap(gsList=None,layerList=None,telDiam=None,telSec=None,fill=False,telCol="red",tells="solid",title=0,outfile=None,sourcedir=None,nx=None,scale=8):
    """
    Shows a plot of guide star overlaps.
    gsList is a list of either LGS objects, NGS objects or tuples of (alt,theta,phi).
    layerList is a list of layer heights
    sourcedir is an optional list of (theta,phi) tuples which shows where, for each layer, the phase is evaluated.
    """
    import pylab
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches
    if gsList==None:
        txt=raw_input("Enter your guide star positions as a list of (alt,theta,phi), e.g. [(90000.,20,x*90) for x in range(4)]: ")
        gsList=eval(txt)
    if layerList==None:
        txt=raw_input("Enter the heights at which to view the overlap (in m) e.g. [0,1000,4000,10000]: ")
        layerList=eval(txt)
    if telDiam==None:
        txt=raw_input("Enter telescope diameter: ")
        telDiam=eval(txt)
        if telSec==None:
            txt=raw_input("Enter telescope secondary diameter: ")
            telSec=eval(txt)
        if fill==False:
            fill=eval(raw_input("Fill?  True/False 1/0: "))
    if telSec!=None:
        telSec/=2.
    
    nlayers=len(layerList)
    if nx==None:
        nx=int(numpy.ceil(numpy.sqrt(nlayers)))
    ny=(nlayers+nx-1)/nx
    f=pylab.figure(figsize=(ny*scale,nx*scale))
    plotno=0
    maxtheta=0
    maxheight=numpy.max(layerList)
    for g in gsList:
        if type(g) in [type(()),type([])]:
            alt,theta,phi=g[:3]
        else:
            alt=g.alt#height
            theta=g.theta
            phi=g.phi
        if theta>maxtheta:
            maxtheta=theta
    lim=maxheight*numpy.tan(maxtheta/3600./180.*numpy.pi)+telDiam/2.
    for l in layerList:
        plotno+=1
        ax=pylab.subplot(nx,ny,plotno)
        ax.add_patch(mpatches.PathPatch(makeDoughnut(telDiam/2.,telSec,0,0),fill=False,ec=telCol,ls=tells))
        #ax.add_patch(pylab.Circle((0,0),telDiam/2.,fill=False,ec="red"))
        #if telSec!=None:
        #    ax.add_patch(pylab.Circle((0,0),telSec/2.,fill=False,ec="red"))
        obs=[]
        for g in gsList:
            if type(g) in [type(()),type([])]:
                alt,theta,phi=g[:3]
                if len(g)>3:col=g[3]
                else:
                    if alt<0:col="blue"
                    elif alt>50000:col="orange"
                    else:col="green"
            else:
                alt=g.alt
                theta=g.theta
                phi=g.phi
                if isinstance(g,LGS):
                    if alt>50000:
                        col="orange"
                    else:
                        col="green"
                elif isinstance(g,NGS):
                    col="blue"
                elif hasattr(g,col):
                    col=g.col
                else:
                    col="blue"
            if alt<0:#ngs
                diam=telDiam
            elif alt<l:#spot formed below layer...
                diam=0
            else:#lgs
                diam=float(alt-l)*telDiam/alt
            r=l*numpy.tan(theta/3600./180.*numpy.pi)
            x=r*numpy.cos(phi/180.*numpy.pi)
            y=r*numpy.sin(phi/180.*numpy.pi)
            if diam>0:
                if telSec==None:
                    sec=None
                else:
                    sec=telSec/telDiam*diam
                ax.add_patch(mpatches.PathPatch(makeDoughnut(diam/2.,sec,x,y),fill=fill,fc=col,alpha=0.25))
                #ax.add_patch(pylab.Circle((x,y),diam/2.,fill=fill,fc=col,alpha=0.25))
                #if telSec!=None:
                #    obs.append((x,y,diam))
        #for o in obs:
        #    ax.add_patch(pylab.Circle((o[0],o[1]),telSec/telDiam*o[2]/2.,fill=fill,fc="white",alpha=0.25))
        pylab.ylim([-lim,lim])
        pylab.xlim([-lim,lim])
        if title:
            pylab.title("%gm"%l)
        if sourcedir!=None:
            for s in sourcedir:
                t,p=s
                r=l*numpy.tan(t/3600./180.*numpy.pi)
                x=r*numpy.cos(p/180.*numpy.pi)
                y=r*numpy.sin(p/180.*numpy.pi)
                ax.add_patch(pylab.Circle((x,y),telDiam/2,fill=False,fc="blue"))
    if outfile!=None:
        pylab.savefig(outfile,bbox_inches="tight")
        pylab.close()
    else:
        pylab.show()

def simple(paramfile=None,batchno=0,lgsAlt=15000.,gateDepth=400.,sig=None,fname=None,lgs_nsubx=None,lgs_clipsize=None,telDiam=None,telSec=None,wfs_int=None,wfs_lat=None,r0=None):
    if fname==None or fname=="None":
        fname="lgsb%da%gd%gs%g.fits"%(batchno,lgsAlt,gateDepth,sig)
    if paramfile!=None:
        config=base.readConfig.AOXml(paramfile,batchno=batchno)
        if lgs_nsubx==None:
            lgs_nsubx=config.getVal("lgs_nsubx")
        if lgs_clipsize==None:
            lgs_clipsize=config.getVal("lgs_clipsize")
        if telDiam==None:
            telDiam=config.getVal("telDiam")
        if telSec==None:
            telSec=config.getVal("telSec")
        if wfs_int==None:
            wfs_int=config.getVal("wfs_int")
        if wfs_lat==None:
            wfs_lat=config.getVal("wfs_lat")
        if r0==None:
            r0=config.getVal("r0")
    npxl=lgs_nsubx*lgs_clipsize
    pup=util.tel.Pupil(npxl,npxl/2.,npxl/2.*telSec/telDiam).fn
    subapFov=5.0
    telFocalLen=46.2
    subapFocalLen=209.
    lgsPower=7.2
    frameRate=1./(wfs_int+wfs_lat)
    gateDepth=gateDepth
    launchAp=0.35
    lgsAlt=lgsAlt
    quality=0.8
    fluxObj=util.guideStar.Flux(wavelength=532.,tel_trans=0.6,optics_trans=1.,launch_trans=1.,pulse_rate=5000.,power=lgsPower,frame_rate=frameRate,QE=1.0)

    lgs=util.guideStar.wfs(wavelength=532.,pup=pup,nsubx=lgs_nsubx,subap_pix=lgs_clipsize,subap_fov=subapFov,LGS_alt=lgsAlt,LGS_gate_depth=gateDepth,launch_ap=launchAp,tel_diam=telDiam,r0=r0,tel_fl=telFocalLen,subap_fl=subapFocalLen,fluxObj=fluxObj,launch_quality=quality,contrast_ratio=5000.)

    lgs.wfs_image()
    mx=max(lgs.spotPsf.sum(3).sum(2).ravel())
    print "Computed photons per subap:",mx
    if sig!=None:
        print "Rescaling to %g photons/subap"%sig
        lgs.spotPsf*=sig/mx
    print "Writing to %s"%fname
    util.FITS.Write(lgs.spotPsf,fname)
    return lgs


class wfsOverview:
    def __init__(self,wfsDict):
        self.wfsDict=wfsDict
    def getWfsByID(self,idstr):
        wfs=self.wfsDict[idstr]
        if wfs.idstr!=idstr:
            raise Exception("Inconsistency in wfsOverview")
        return wfs
    def values(self):
        return self.wfsDict.values()


"""
getWfsByID(idstr)


wfsobj needs:
nsubx
nimg
phasesize (wfs_n)
minarea
pupil
seed (None)
#calSource (0)
integSteps  (wfs_int/tstep)
atmosPhaseType (phaseonly)   ["phaseonly","phaseamp","realimag"]
nfft  (phasesize*2)
clipsize (nfft)
nimg  (clipsize/2)
bglevel (wfs_read_mean)
readoutNoise (wfs_read_sigma)
rowint (None)
threshType (0)
latency
skyBrightness
ncen (nimg)
floor
sig
addPoisson
lgsPsf (was laserGuideStar)
spotpsf
opticalBinning
magicCentroiding
linearSteps
calNCoeff
stepRangeFrac
centWeight (None)
correlationCentroiding (0)
corrThresh
corrPattern (None)
useBrightest (0)
preBinningFactor (1)
sourcelam
"""

class LGS(util.atmos.source):
    """Simple structure for holding LGS info.
    Important parameters are:
    nsubx - number of subaps
    theta, phi - radial coords of on-sky positions, in arcsec, degrees.
    height - the height of the laser emission (height of guide star).
    phasesize - number of phase pixels in 1 sub-aperture.
    minarea - fraction of vignetting before sub-aperture is no longer used.
    sig - WFS signal in photons per subap per frame.
    idstr - the identification string of the light path of this wfs.
    sourcelam - wavelength (nm) of the guide star.
    phslam - wavelength (nm) of the wavefront phase.
    reconList - list of reconstructors which use this slope information.
    pupil - the telescope pupil function for this guide star.
    fov - for widefield systems - radius, not diam
    """
    def __init__(self,idstr,nsubx,theta,phi,height,phasesize,pupil,minarea=0.5,sig=1e6,launchDist=0.,launchTheta=0.,sourcelam=None,phslam=None,reconList=None,nimg=None,nfft=None,clipsize=None,ncen=None,preBinningFactor=1,bglevel=0.,readoutNoise=0.,integSteps=1,rowint=None,threshType=0,latency=0,skyBrightness=0.,floor=0.,seed=0,atmosPhaseType="phaseonly",addPoisson=1,lgsPsf=None,spotpsf=None,opticalBinning=0,magicCentroiding=0,linearSteps=None,calNCoeff=0,stepRangeFrac=1,centWeight=None,correlationCentroiding=0,corrThresh=0,corrPattern=None,useBrightest=0,fov=0.,parabolicFit=0,gaussianFitVals=None,subapFlag=None,integstepFn=None):
        """pupil can be a util.tel object"""
        #Initialise the parent...
        super(LGS,self).__init__(idstr,theta,phi,height,sourcelam=sourcelam,phslam=phslam,sig=sig)
        #self.theta=theta#angle in arcsec
        #self.phi=phi#angle in degrees
        self.nsubx=nsubx
        self.phasesize=phasesize
        self.nimg=nimg
        self.nfft=nfft
        self.clipsize=clipsize
        self.ncen=ncen
        self.fov=fov#radius not diam
        self.preBinningFactor=preBinningFactor
        self.bglevel=bglevel#was wfs_read_mean
        self.readoutNoise=readoutNoise#was wfs_read_sigma
        self.integSteps=integSteps#in units of tstep
        self.rowint=rowint#in units of tstep
        self.threshType=threshType
        self.latency=latency#in units of tstep
        self.skyBrightness=skyBrightness
        self.floor=floor
        self.seed=seed
        self.atmosPhaseType=atmosPhaseType
        self.addPoisson=addPoisson
        self.lgsPsf=lgsPsf
        self.spotpsf=spotpsf
        self.opticalBinning=opticalBinning
        self.magicCentroiding=magicCentroiding
        self.linearSteps=linearSteps
        self.calNCoeff=calNCoeff
        self.stepRangeFrac=stepRangeFrac
        self.centWeight=centWeight
        self.correlationCentroiding=correlationCentroiding
        self.corrThresh=corrThresh
        self.corrPattern=corrPattern
        self.parabolicFit=parabolicFit
        self.gaussianFitVals=gaussianFitVals#None, or a tuple of (gaussianMinVal, gaussianReplaceVal)
        self.useBrightest=useBrightest
        #self.alt=height Now self.alt.  
        self.minarea=minarea
        #self.sig=sig
        #self.idstr=idstr
        self.reconList=reconList
        self.launchDist=launchDist#distance of laser launch off-axis, in m.
        self.launchTheta=launchTheta#angle of launch from 3 oclock, anticlockwise
        self.integstepFn=integstepFn
        self.dmheight=None#used in computeCoords, to avoid recomputation if not necessary.
        self.coords=None#used in computeCoords
        self.telDiam=None#used in computeCoords
        self.pupil=pupil
        if type(pupil)==numpy.ndarray:
            pup=pupil
        else:
            pup=pupil.fn
        if pup.shape[0]!=nsubx*phasesize:
            raise Exception("Wrong sized pupil: %d (should be %d)"%(pup.shape[0],nsubx*phasesize))
        self.subapFlag=subapFlag#use getSubapFlag() to get this.
        if self.phasesize==None:
            if self.nimg!=None:
                self.phasesize=self.nimg
        if self.phasesize!=None:
            if self.nfft==None:
               self.nfft=self.phasesize*2
            if self.clipsize==None:
                self.clipsize=self.nfft
            if self.nimg==None:
                self.nimg=self.clipsize//2
            if self.ncen==None:
                self.ncen=self.nimg

    def computeCoords(self,telDiam,height,fname=None):
        """Creates an array containing coords of centre of each subap.
        Height is the height of interest (eg the DM conjugate height), not the height of the LGS.
        0,0 means on axis."""
        if self.coords!=None and telDiam==self.telDiam and height==self.dmheight:
            return#already computed
        self.telDiam=telDiam
        if self.alt<=height:
            raise Exception("LGS height lower than DM height %g %g"%(self.alt,height))
        telDiam*=(self.alt-height)/float(self.alt)
        self.dmheight=height
        self.coords=numpy.zeros((self.nsubx,self.nsubx,2),numpy.float32)
        r=height*numpy.tan(self.theta/3600./180.*numpy.pi)
        wsx=r*numpy.cos(self.phi/180.*numpy.pi)#position of centre of phase on dm.
        wsy=r*numpy.sin(self.phi/180.*numpy.pi)
        #allow for cone effect
        subapDiam=telDiam/self.nsubx#*(self.alt-height)/float(self.alt)
        if fname!=None:
            f=open(fname,"a")
            f.write("#GS height %g, nsubx %g, theta %g, phi %g\n"%(height,self.nsubx,self.theta,self.phi))
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                self.coords[i,j,0]=wsx-telDiam/2.+j*subapDiam+subapDiam/2
                self.coords[i,j,1]=wsy-telDiam/2.+i*subapDiam+subapDiam/2
                if fname!=None:
                    f.write("%g\t%g\n"%(self.coords[i,j,0],self.coords[i,j,1]))
        if fname!=None:
            f.write("\n")
            f.close()
    def getSubapFlag(self):
        """Compute the subap flags for a given nsubx"""
        if self.subapFlag!=None:
            return self.subapFlag
        nsubx=self.nsubx
        subflag=numpy.zeros((nsubx,nsubx),numpy.int32)
        n=self.phasesize#self.npup/nsubx
        minarea=self.minarea*n*n
        if type(self.pupil)==numpy.ndarray:
            pup=self.pupil
        else:
            pup=self.pupil.fn
        for i in range(nsubx):
            for j in range(nsubx):
                subarea=pup[i*n:(i+1)*n,j*n:(j+1)*n].sum()
                if subarea>=minarea:
                    subflag[i,j]=1
        return subflag


class NGS(util.atmos.source):
    """simple structure for holding NGS info
    Important parameters are:
    nsubx - number of subaps
    theta, phi - radial coords of on-sky positions, in arcsec, degrees.
    phasesize - number of phase pixels in 1 sub-aperture.
    pupil - the pupil function.
    minarea - fraction of vignetting before sub-aperture is no longer used.
    sig - WFS signal in photons per subap per frame.
    idstr - the identification string of the light path of this wfs.
    sourcelam - wavelength (nm) of the guide star.
    phslam - wavelength (nm) of the wavefront phase.
    reconList - list of reconstructors which use this slope information.
    pupil - the telescope pupil function for this guide star.
    fov - radius of the fov, not diam, in arcsec
    """
    def __init__(self,idstr,nsubx,theta,phi,phasesize,pupil,minarea=0.5,sig=None,sourcelam=None,phslam=None,reconList=None,nimg=None,nfft=None,clipsize=None,ncen=None,preBinningFactor=1,bglevel=0.,readoutNoise=0.,integSteps=1,rowint=None,threshType=0,latency=0,skyBrightness=0,floor=0.,seed=0,atmosPhaseType="phaseonly",addPoisson=1,spotpsf=None,opticalBinning=0,magicCentroiding=0,linearSteps=None,calNCoeff=0,stepRangeFrac=1,centWeight=None,correlationCentroiding=0,corrThresh=0,corrPattern=None,useBrightest=0,fov=0.,parabolicFit=0,gaussianFitVals=None,subapFlag=None,integStepFn=None):
        """pupil can be a util.tel object"""
        #Initialise parent...
        super(NGS,self).__init__(idstr,theta,phi,-1,sourcelam=sourcelam,phslam=phslam,sig=sig)

        #self.theta=theta#angle in arcsec
        #self.phi=phi#angle in degrees
        self.nsubx=nsubx
        self.phasesize=phasesize
        self.nimg=nimg
        self.nfft=nfft
        self.clipsize=clipsize
        self.ncen=ncen
        self.fov=fov#radius, not diam.
        self.preBinningFactor=preBinningFactor
        self.bglevel=bglevel#was wfs_read_mean
        self.readoutNoise=readoutNoise#was wfs_read_sigma
        self.integSteps=integSteps#in units of tstep
        self.rowint=rowint#in units of tstep
        self.threshType=threshType
        self.latency=latency#in units of tstep
        self.skyBrightness=skyBrightness
        self.floor=floor
        self.seed=seed
        self.atmosPhaseType=atmosPhaseType
        self.addPoisson=addPoisson
        self.lgsPsf=None
        self.spotpsf=spotpsf
        self.opticalBinning=opticalBinning
        self.magicCentroiding=magicCentroiding
        self.linearSteps=linearSteps
        self.calNCoeff=calNCoeff
        self.stepRangeFrac=stepRangeFrac
        self.centWeight=centWeight
        self.correlationCentroiding=correlationCentroiding
        self.corrThresh=corrThresh
        self.corrPattern=corrPattern
        self.parabolicFit=parabolicFit
        self.gaussianFitVals=gaussianFitVals#None, or a tuple of (gaussianMinVal, gaussianReplaceVal)
        self.useBrightest=useBrightest
        self.minarea=minarea
        #self.sig=sig
        #self.idstr=idstr
        self.integstepFn=integStepFn
        self.reconList=reconList
        self.dmheight=None#used in computeCoords
        self.coords=None#used in computeCoords
        self.telDiam=None#used in computeCoords
        self.pupil=pupil
        self.subapFlag=subapFlag#use getSubapFlag() to get this.
        if self.phasesize==None:
            if self.nimg!=None:
                self.phasesize=self.nimg
        if self.phasesize!=None:
            if self.nfft==None:
                self.nfft=self.phasesize*2
            if self.clipsize==None:
                self.clipsize=self.nfft
            if self.nimg==None:
                self.nimg=self.clipsize//2
            if self.ncen==None:
                self.ncen=self.nimg

    def computeCoords(self,telDiam,height,fname=None):
        """Creates an array containing coords of centre of each subap.
        Height is the height of interest (eg the DM conjugate height).
        0,0 means on axis."""
        if self.coords!=None and telDiam==self.telDiam and height==self.dmheight:
            return#already computed
        self.telDiam=telDiam
        self.dmheight=height
        self.coords=numpy.zeros((self.nsubx,self.nsubx,2),numpy.float32)
        r=height*numpy.tan(self.theta/3600./180.*numpy.pi)
        wsx=r*numpy.cos(self.phi/180.*numpy.pi)#position of centre of WFS (0,0=onaxis) in m
        wsy=r*numpy.sin(self.phi/180.*numpy.pi)
        subapDiam=telDiam/self.nsubx
        if fname!=None:
            f=open(fname,"a")
            f.write("#GS height %g, nsubx %g, theta %g, phi %g\n"%(height,self.nsubx,self.theta,self.phi))
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                self.coords[i,j,0]=wsx-telDiam/2.+j*subapDiam+subapDiam/2
                self.coords[i,j,1]=wsy-telDiam/2.+i*subapDiam+subapDiam/2
                if fname!=None:
                    f.write("%g\t%g\n"%(self.coords[i,j,0],self.coords[i,j,1]))
        if fname!=None:
            f.write("\n")
            f.close()

    def getSubapFlag(self):
        """Compute the subap flags for a given nsubx"""
        if self.subapFlag!=None:
            return self.subapFlag
        nsubx=self.nsubx
        subflag=numpy.zeros((nsubx,nsubx),numpy.int32)
        n=self.phasesize#self.npup/nsubx
        minarea=self.minarea*n*n
        if type(self.pupil)==numpy.ndarray:
            pup=self.pupil
        else:
            pup=self.pupil.fn
        for i in range(nsubx):
            for j in range(nsubx):
                subarea=pup[i*n:(i+1)*n,j*n:(j+1)*n].sum()
                if subarea>=minarea:
                    subflag[i,j]=1
        return subflag


#Version 3 - analysis
#Includes contrast ratio - Not Tested
#Includes cross-talk between subapertures - Not Tested
#Includes residual spot offset calculations - Tested
#Accurate LGS analysis not implemented


class wfs:
    """A class to compute LGS spot elongation.  Call wfs_image() to get the spots.
    Originally from TJM.  Needs tidying up and speeding up.
    Notes for ELT:
    Need to increase LGS_vert_res to eg 30.
    tel_diam and tel_fl need adjusting
    nsubx needs adjusting
    subap_fov probably needs adjusting
    subap_pix should probably also be changed?
    Fairly nice spots can be obtained using:
    wfs(tel_diam=42.,tel_fl=200.,nsubx=32,subap_fov=3.,r0=0.2,subap_pix=16,LGS_vert_res=10.,LGS_gate_depth=500.,launch_ap=0.5)

    For sodium, could try:
    wfs=util.guideStar.wfs(LGS_alt=92000.,nsubx=32,subap_fov=14.,tel_diam=42.,wavelength=589.,LGS_gate_depth=25000.,LGS_vert_res=20,r0=0.2,tel_fl=600.,subap_pix=16,launch_ap=0.5,fluxObj=util.guideStar.Flux(power=5.,frame_rate=800.))
    Note, subap_fov = nsubx * lenslet_pitch**2 / ( 2 * tel_diam * lenslet_fl )
    """

    def __init__(self,pup=None,nsubx=8,subap_pix=8,subap_fov=2.8,tel_diam=4.2,wavelength=None,LGS_alt=34000.,max_vignette=0.3,LGS_gate_depth=5000.,LGS_vert_res=10.,launch_quality=0.8,r0=0.11,launch_ap=0.5,trunc=0.83,m_squared=1.2,tel_fl=46.2,subap_pitch=None,subap_fl=209.,LGS_theta=0.,LGS_phi=0.,gauss=1,tel_sec=None,wfs_addnoise=0,wfs_read=10.,wfs_floor=None,contrast_ratio=500.,fluxObj=None,lgsDefault=None,laserType=None):
        """pupfn is the pupil transmission function (util.tel.Pupil)
        nsubx is number of subaps in 1 direction
        subap_pix is number of pixels to use per subap (wfs_clipsize)
        subap_fov is field of view
        tel_diam is telescope primary diameter
        wavelength is the wavelength to use, nm
        LGS_alt is LGS altitude
        max_vignette is minimum lenslet illumination ratio
        LGS_gate_depth is height of LGS range rate in m
        LGS_vert_res is number of bins to divide this into.
        launch_quality is ???
        r0 is Frieds parameter
        launch_ap is ??? launch telescope aperture size???
        trunc is ???
        m_squared is laser quality measurement
        tel_fl is telescope focal length in m
        subap_pitch is pitch of subaps in mm
        subap_fl is in mm
        LGS_theta, LGS_phi are sky position of LGS.
        gauss is flag, whether to use gaussian or sinc LGS profile.
        tel_sec is telescope secondary diameter, used if pupfn not specified.
        wfs_addnoise is whether to add noise
        wfs_read is the readout noise for WFS.
        wfs_floor is the floor to subtract from WFS
        contrast_ratio is the shutter contrast ratio.
        fluxObj is a Flux object.
        lgsDefault if not None, can be used to automatically select eg wavelength, alt etc.
        """
        defaults={"sodium":{"wavelength":589.,"LGS_alt":92000.,"LGS_gate_depth":25000.,"LGS_vert_res":20.,},
                  "rayleigh":{"wavelength":532.,"LGS_alt":34000.,"LGS_gate_depth":5000.,"LGS_vert_res":10.,},
                  }
        if lgsDefault in defaults.keys():
            wavelength=defaults[lgsDefault]["wavelength"]
            LGS_alt=defaults[lgsDefault]["LGS_alt"]
            LGS_gate_depth=defaults[lgsDefault]["LGS_gate_depth"]
            LGS_vert_res=defaults[lgsDefault]["LGS_vert_res"]
            print "LGS_alt %g, gate depth %g wavelength %g, resolution %d"%(LGS_alt,LGS_gate_depth,wavelength,LGS_vert_res)
        self.spotPsf=None
        self.subapImage=None#will hold the elongated spot array eventually
        self.sig=0.#will hold the photons in each subap eventually.
        #paramfile ='params.py'
        self.tel_diam=tel_diam
        #self.tel_alt=tel_alt
        self.subap_pix=subap_pix
        self.subap_fov=subap_fov
        self.nsubx=nsubx
        self.LGS_alt=LGS_alt
        if laserType==None:
            if LGS_alt>75000.:
                print "INFORMATION Assuming Sodium laser"
                self.laserType="sodium"
            else:
                print "INFORMATION Assuming Rayleigh laser"
                self.laserType="rayleigh"
        else:
            self.laserType=laserType
        npup=subap_pix*nsubx
        if type(pup)==type(None):
            if tel_sec==None:
                tel_sec=tel_diam/7.
            pup=util.tel.Pupil(npup,npup/2,npup/2.*tel_sec/tel_diam).fn
        if pup.shape!=(npup,npup):
            if type(pup)==numpy.ndarray:
                raise Exception("Pupil function for util.guidestar is wrong shape (%s), should be %s"%(str(pup.shape),str((npup,npup))))
            else:#create a scaled up pupil...
                print "INFORMATION util.guideStar - scaling up pupil function for LGS spots"
                pup=util.tel.Pupil(npup,pup.r1*npup/pup.npup,pup.r2*npup/pup.npup).fn
        self.pup=pup
        self.max_vignette=max_vignette
        self.LGS_gate_depth=LGS_gate_depth
        self.LGS_vert_res=LGS_vert_res
        self.launch_quality=launch_quality
        self.r0=r0
        if wavelength==None:
            if wavelength=="sodium":
                wavelength=589.
            elif wavelength=="rayleigh":
                wavelength=514.
        self.wavelength=wavelength
        self.launch_ap=launch_ap
        self.trunc=trunc
        self.m_squared=m_squared
        self.tel_fl=tel_fl
        if subap_pitch==None:#compute subap pitch...
            subap_pitch=numpy.sqrt(subap_fov/3600./180*numpy.pi * 2 * tel_diam * subap_fl/1000. / nsubx)*1000.
            print "subap_pitch now %g mm"%subap_pitch
        if subap_fl==None:#compute subap focal length
            subap_fl=nsubx*(subap_pitch/1000.)**2/(2*tel_diam*subap_fov/3600./180*numpy.pi)*1000.
            print "subap_fl now %g mm"%subap_fl
        if subap_fov==None:#compute subap fov
            subap_fov=nsubx*(subap_pitch/1000.)**2/(2*tel_diam*(subap_fl/1000.))*180*3600/numpy.pi
            
        #Now check that the pitch, fov and focal lengths are consistent...
        realfov=nsubx*(subap_pitch/1000.)**2/(2*tel_diam*(subap_fl/1000.))*180*3600/numpy.pi
        if numpy.fabs(realfov-subap_fov)>1e-6:
            raise Exception("Disagreement between subap pitch, focal length and fov (fov specified as %g, calculated as %g"%(subap_fov,realfov))
        self.subap_pitch=subap_pitch
        self.subap_fl=subap_fl
        self.LGS_theta=LGS_theta
        self.LGS_phi=LGS_phi
        self.gauss=gauss
        self.wfs_addnoise=wfs_addnoise
        self.wfs_read=wfs_read
        if wfs_floor==None:
            wfs_floor=wfs_read*3
        self.wfs_floor=wfs_floor
        self.contrast_ratio=contrast_ratio
        self.wfs_scale = self.subap_fov/float(self.subap_pix)
        #Generates blank output wfs image
        self.npup = self.nsubx*self.subap_pix
        self.pxlArea=(self.tel_diam/float(self.npup))**2
        #print npup
        self.wfs = numpy.zeros((self.npup,self.npup)).astype('f')
        # Generate WFS pupil mask
        #obs = (self.tel_2diam/float(self.tel_diam))*self.npup/2.
        #self.pup = self.pupil(self.npup,self.npup/2.,obs)
        
        # Generate positions of lenslets in terms of WFS pixels
        subaps = self.subap_pos(self.nsubx,self.subap_pix)

        # Calculate lenslet properties - square lenslets
        lenslet = []
        n = int((self.subap_pix/2.)+0.5)
        min_dist = 1E31
        #closest_lenslet = 0
        for i in range(subaps.shape[0]):
            x = subaps[i,0]
            y = subaps[i,1]
            # Need to find closest lenslet to determine plume distance to calculate
            #if ((x**2)+(y**2))<min_dist: closest_lenslet=i
            # Now need to find out light collecting area of each lenslet
            pup_slice = self.pup[x-n:x+n,y-n:y+n]
            a = numpy.sum(numpy.sum(pup_slice))*self.pxlArea
            lenslet.append((x,y,a))
        # self.lenslet contains lenslet x and y positions and collecting area in m^2
        self.lenslet=numpy.array(lenslet)
        #Also need to know maximum filled lenslet area
        max_area = (self.subap_pix**2)*self.pxlArea
        self.subap_calc_limit = max_area * self.max_vignette

        # Initialise Rayleigh photon return - required even for sodium LGS
        if fluxObj==None:
            print "WARNING: no flux object given for LGS - assuming default parameters"
            self.flux = Flux(wavelength=wavelength)#wavelength,tel_trans,optics_trans,power,launch_trans,pulse_rate,frame_rate,QE,zenith,tel_alt,scale_height)
        else:
            self.flux=fluxObj

        # Calculate maximum and minimum plume heights
        # Use closest lenslet as this has the maximum depth of field

	#offset = sqrt(self.lenslet[closest_lenslet][0]**2+self.lenslet[closest_lenslet][1]**2) #pixels
        #offset = offset*self.tel_diam/float(self.npup) #metres
        #max_dist = 1.5*self.subap_fov*sqrt(2) #1/sqrt(2) = cos(45)
        ## max dist assumes that the spot is elongated at 45 degrees and the plume needs
        ## to be calculated across a subaperture 2xsubap_pix across
        ## This value is the absolute maximum
        #column_upper,column_lower = self.plume_depth(offset,max_dist,self.LGS_alt)

        column_upper = (self.LGS_gate_depth/2.)+self.LGS_alt
        column_lower = self.LGS_alt-(self.LGS_gate_depth/2.)

        # We know we want self.LGS_vert_res steps across range gate to avoid z-granularity
        # Outside range gate depth, we still want to have same vertical sampling
        self.LGS_steps = int((((column_upper-column_lower)/float(self.LGS_gate_depth))*self.LGS_vert_res)+0.5)
        #print column_upper, column_lower, self.LGS_gate_depth, self.LGS_vert_res
        print "Using %g LGS_steps"%self.LGS_steps
        #raw_input('N:')
        depth = (column_upper-column_lower)/float(self.LGS_steps)
        # We now know depth of each step through the plume and how many steps through entire plume
        # print depth, self.LGS_steps, column_lower, column_upper, column_lower+((self.LGS_steps+0.5)*depth)
        # Now calculate LGS properties for each step
        rgd = self.LGS_gate_depth
        self.lgs_return_data = numpy.zeros((self.LGS_steps,3),numpy.float64)
        if self.laserType=="sodium":
            sodiumReturn=self.flux.sodium(range(int(column_lower+0.5*depth),int(column_lower+(self.LGS_steps+0.5)*depth),int(depth)))
        for i in range(self.LGS_steps):
            alt = column_lower+((i+0.5)*depth)
            self.lgs_return_data[i,0] = alt
            if self.laserType=="sodium":
                self.lgs_return_data[i,1]=sodiumReturn[i]
            else:
                self.lgs_return_data[i,1] =self.flux.rayleigh(alt,depth)#self.lgs_return(alt,depth) #This is in photons per m^2
            self.lgs_return_data[i,2] = self.lgs_diam(alt)
            test = numpy.fabs(self.LGS_alt-alt)
            if test>=(rgd/2.):
                self.lgs_return_data[i,1]=self.lgs_return_data[i,1]/float(self.contrast_ratio)
            #else:
                #print alt, 'is inside range gate'

    def plume_depth(self,offs,dist,alt):
        lgs_ang = numpy.arctan(dist/float(alt))*206265
        min_ang = lgs_ang+offs
        max_ang = lgs_ang-offs
        min_alt = dist/tan(min_ang/206265.)
        max_alt = dist/tan(max_ang/206265.)
        return max_alt, min_alt

    def blt_res(self,atmosBlurring=1):
        """calculate fwhm of lgs spot"""
        # First calculate radians^2 of error in optical path
        var_qual = -numpy.log(self.launch_quality)
        # Calculate r0 at laser wavelength
        r0_las = self.r0*((self.wavelength/500.)**(6/5.))
        if atmosBlurring:
            # Now calculate short-exposure Strehl ratio
            rho = float(r0_las*(1+(0.37*(r0_las/self.launch_ap)**(1/3.))))
        else:
            rho=1.
        #ex = exp(-qual)
        Drho = float(1.+(self.launch_ap/rho)**2)
        # Calculate variance due to tilt removed turbulence
        var_tilt = 0.134*((self.launch_ap/float(r0_las))**(5/3.))
        # Now we have total variance of beam
        var_tot = var_tilt+var_qual
        # Add factor for gaussian beam propagation
        K = 1.6449+(0.6460/(self.trunc-0.2816)**1.821)-(0.5320/(self.trunc-0.2816)**1.891)
        # Now calculate short-exposure resolution
        strehl = numpy.exp(-var_tot)+(1-numpy.exp(-var_tot))/Drho
        #print r0_las,rho,var_tot,Drho,strehl,K
        #res = K*(1.22*self.wavelength*1E-9/self.launch_ap)*((qual*qual)+((1-qual)**2/Drho))**0.5
        res = K*(1.22*self.wavelength*1E-9/self.launch_ap/strehl)*(numpy.exp(var_tot*-2)+((1-numpy.exp(-var_tot))**2/Drho))**0.5
        # Add effect of M-squared
        res = res * numpy.sqrt(self.m_squared)*self.LGS_alt
        #res = 0.048481/2. #This is 1 arcsecond diameter (0.5" radius) LGS at 20km
        #angle = self.LGS_min_diam/2./206265. #arcseconds to radians
        #res = self.LGS_alt*tan(angle)
        #print 'res:',res
        return res

    '''def tel_focus(self,alt_focus,alt_slice,tel_fl):
        # Calculate focal lengths of beams coming from different altitudes inside telescope
        LGS_fl = 1/float((1/float(tel_fl))-(1/float(alt_focus)))
        slice_fl = 1/float((1/float(tel_fl))-(1/float(alt_slice)))
        # Now we have focal point difference between WFS focal plane and focal point of slice
        # Want to determine relative magnifications of the two beams
        # Can compare focal lengths of beams on-sky and in telescope to get magnification
        #magnification_slice = alt_slice/slice_fl      # Magnification of beam coming from slice altitude
        #magnification_LGS = alt_focus/LGS_fl          # Magnification of beam coming from LGS altitude
        #ratio_slice = magnification_slice/magnification_LGS
        # Now we need to determine actual increase in diameter due to defocus on WFS
        # There are lenslet arrays etc, so we just want to know the ratio with which to
        # increase the LGS diameter by, as this is independent of post-telescope optics
        tel_fr_slice = slice_fl/float(self.tel_diam)
        diameter = float(abs(LGS_fl-slice_fl))/float(tel_fr_slice)
        # This is the diameter of the defocused image in meters
        # alt_slice,tel_fl and tel_diam should all be in the same units

        # 1 arcsecond at slice_alt meters from telescope aperture is
        arcsec = alt_slice/206265. # This is in meters
        # Telescope magnification is given by slice_fl/alt_slice (thin lens approximation)
        mag_slice = slice_fl/alt_slice
        # Therefore the plate scale of the telescope when looking at a point slice_alt meters away is
        scale_slice = arcsec*mag_slice #This is in meters
        demag = diameter/scale_slice #This is angle equal to the defocus diameter in arcseconds
        res = alt_slice*tan(demag/2./206265.) #Now have radius of defocus in metres at slice altitude
        # This is added to the blt_res to simulate the effect of the single conjugate plane of the WFS CCD
        # The diameter of the LGS plume at the slice altitude is effectively increased by a factor res
        #return LGS_fl,slice_fl,magnification_slice,magnification_LGS,ratio_slice,tel_fr_slice,diameter,ratio_LGS,scale_LGS,scale_slice,demag
        return res'''

    def tel_focus(self,alt_focus,alt_slice,tel_fl,tel_diam,lenslet_pitch,n_lenslets,lenslet_fl):
        """
        Calculates defocus of subapertures on WFS CCD plane
        using a 3 lens (telescope,collimator,lenslet) thin-lens approximation to the WFS train
        Input parameters are:
            alt_focus     = the focal altitude of the LGS
            alt_slice     = the focal altitude of the slice in the Rayleigh plume
            tel_fl        = telescope focal length
            lenslet_pitch = self explanatory
            n_lenslets    = number of lenslets across aperture
        """
        # First need to calculate a set of constants that are set by the WFS being conjugate to alt_focus
        L2 = 1/((1/float(tel_fl))-(1/float(alt_focus)))  # Plane in telescope primary conjugate focus to altitude, alt_focus
        F2 = n_lenslets*L2*lenslet_pitch/float(tel_diam) # Focal length of collimator to give collimated beam lenslets*pitch in diameter in mm
        F23 = 1/((1/F2)+(1/float(lenslet_fl)))           # Focal length of collimator and lenslet combined
        LT = L2+F2
        #print L2,F2,F23,LT
        # Now we can calculate the actual focal points of the plume at alt_slice
        L3 = float(LT - 1./((1/float(tel_fl))-(1/float(alt_slice)))) # This gives us the distnce between focal points of alt_focus and alt_slice
        # Put it through the collimator and lenslet
        L4 = 1/((1/F23)-(1/L3))
        # Calculate f-ratio of output from lenslet
        output_fr = L4/float(lenslet_pitch)
        #print L3,L4,output_fr
        # Now can caculate actual defocused image diameter - the increase in diameter.
        diameter = (L4-lenslet_fl)/output_fr
        # Can work out magnification between input and output beams
        input_fr = alt_slice*n_lenslets/float(tel_diam)  # We want this in terms of a single lenslet so divide the telescope aperture by n_lenslets
        # System magnification is given by comparing input and output f-ratios
        mag = input_fr/output_fr
        #print diameter,input_fr,mag
        # Therefore the diameter is magnified by this amount on sky
        res = diameter*mag
        #Now return the (magnitude of) the increase in arcseconds of the size of the slice diameter
        return abs(res)

    '''def tel_focus3(self,alt_focus,alt_slice,tel_fl,tel_diam,lenslet_pitch,n_lenslets,lenslet_fl):
        # Calculates defocus of subapertures on WFS CCD plane
        # using a 3 lens (telescope,collimator,lenslet) thin-lens approximation to the WFS train
        # Input parameters are:
        #    alt_focus     = the focal altitude of the LGS
        #    alt_slice     = the focal altitude of the slice in the Rayleigh plume
        #    tel_fl        = telescope focal length
        #    lenslet_pitch = self explanatory
        #    n_lenslets    = number of lenslets across aperture
     
        # First need to calculate a set of constants that are set by the WFS being conjugate to alt_focus
        L2 = 1/((1/float(tel_fl))-(1/float(alt_focus)))  # Plane in telescope primary conjugate focus to altitude, alt_focus
        F2 = n_lenslets*L2*lenslet_pitch/float(tel_diam) # Focal length of collimator to give collimated beam lenslets*pitch in diameter
        LT = L2+F2
        #Need to place lenslet conjugate pupil plane
        LC = 1/((1/float(F2))-(1/float(LT)))
        
        print L2,F2,LC,LT
        
        A = array(((1,0),(1/float(lenslet_fl),1)))
        B = array(((1,-LC),(0,1)))
        C = array(((1,0),(1/float(F2),1)))
        D = array(((1,-LT),(0,1)))
        E = array(((1,0),(1/float(tel_fl),1)))
        # Now we can calculate the actual focal points of the plume at alt_slice'''

    def lgs_diam(self,alt):
        # calculate fwhm of lgs spot
        w0 = self.blt_res()
        tmp=w0
        # add defocus term from tel_focus (makes the spot a bit bigger because its not focussed at alt, rather is focussed at lgs_alt.
        w0 = w0 + self.tel_focus(self.LGS_alt,alt,self.tel_fl,self.tel_diam,self.subap_pitch,self.nsubx,self.subap_fl)
        #print w0
        # want an equation that has the form of x = A[B+(Cy)**2]**0.5
        # (identical to Melles-Griot gaussian beam waist optics propagation eqn)
        # that returns a waist of the BLT diameter at 0km and LGS spot size at LGS focal height
        if w0>self.launch_ap:
            print self.launch_ap,w0,tmp
        C = numpy.sqrt(((self.launch_ap/w0)**2.)-1.)
        y = ((self.LGS_alt-alt)/float(self.LGS_alt)) # strictly (self.LGS_alt - alt) always must be positive, but as y is squared this isn't an issue
        B = 1.
        A = w0
        x = A*(B+((C*y)**2))**0.5 #This is the beam radius returned in metres, need a diameter in arcseconds
        #waist = atan(x/float(self.LGS_alt))*206265*2
        waist = numpy.arctan(x/float(alt))*206265*2
        # Use self.LGS_alt unless focal length of lenslet is known, then use alt
        # If focal length of lenslet is not known, the above equation emulates the defocus present in the lenslet
        # by producing a beam waist, otherwise the dominant factor in spot size on the CCD is the angular size of the spot diminishing as y increases
        #print 'waist:',waist, alt
        # Current version of gauss expects fwhm
        # by multiplying by 1.17 and dividing by 2, we get 1/e^2 diameter
        # this is compatible with GLAS BLT aperture predictions calculated by ron humphreys
        # Long exposure approximation for speed
        #arcsecs = 1.
        return waist*1.17/(2.*float(self.wfs_scale))


    def subap_pos(self,num_subap, pixels_per_ap):
        """Generate positions of lenslets in terms of WFS pixels
        """
        ap_posn = numpy.zeros((num_subap**2,2),numpy.int32)
        half_ap = int((pixels_per_ap/2.)+0.5)
        x = half_ap
        y = half_ap
        for i in range(num_subap**2):
            x += pixels_per_ap
            if i%num_subap == 0:
                y += pixels_per_ap
                x = half_ap
            if i == 0:
                y = half_ap
                x = half_ap
            ap_posn[i,:,] = (y,x)
        return ap_posn.astype(numpy.int32)


##############################################################################################

    def ang_offset(self,alt,lgs_alt,offset):
        ang1 = numpy.arctan(lgs_alt/offset)*206265.
        ang2 = numpy.arctan(alt/offset)*206265.
        ang = ang1-ang2
        return ang


    def direction(self,x,y):
        if x>0:
            if y>=0:
                return numpy.arctan(y/x)
            else:
                return numpy.arctan(y/x)+(2*numpy.pi)
        elif x==0:
            if y>0:
                return numpy.pi/2.
            elif y==0:
                return 0
            else:
                return 1.5*numpy.pi
        else:
            return numpy.arctan(y/x)+numpy.pi

##     def mask_image(self,image):
##         # Makes guard ring around subaperture
##         out = numpy.zeros((image.shape),float64)
##         out[1:-2,1:-2]=image[1:-2,1:-2]
##         #image[:,0] = 0.
##         #image[:,-1] = 0.
##         #image[0,:,] = 0.
##         #image[-1,:,] = 0.
##         return out

    def add_noise(self,image):
        if not self.wfs_addnoise:
            return image
        read = self.wfs_read
      	bimg = image.astype(numpy.float64)
        totsig = numpy.sum(numpy.sum(bimg))	
        if(totsig>0.):
            #img0=((bimg/totsig)*self.sig)+read*read
            print "ERROR IN GUIDESTAR: DODGY NOISE ADDITION!!!"
            img0=bimg+read*read
            cmod.imgnoise.shot(img0,bimg)#inject shot noise

            bimg=clip(bimg,self.wfs_floor,1.e6)
            bimg = bimg-self.wfs_floor
        return bimg

    def centroid(self,image):
        im = image.astype(numpy.float64)
        if numpy.sum(numpy.sum(im))==0:
            x=y=4.5
        else:
            a =numpy.sum(im)
            b =numpy.sum(im,1)
            a1 = a * (arange(len(a))+1)
            b1 = b * (arange(len(b))+1)
            x = numpy.sum(a1)/float(numpy.sum(a))
            y = numpy.sum(b1)/float(numpy.sum(b))
        return x,y

    def subap_im(self,x,y,data,pixels,out):
        diam = data[2]
        flux = data[1]
        if flux>0.0:
            # Generate a gaussian centred at x,y with a fwhm of diam scaled to flux
            if self.gauss==1:
                self.offset_gauss(pixels,diam,x,y,out)
            else:
                self.offset_sinc(pixels,diam,x,y,out)
            #Rescale to unity
            #out = out-numpy.minimum.reduce(numpy.minimum.reduce(out))
            #out = out/(numpy.maximum.reduce(numpy.maximum.reduce(out))+1E-10)
            #Find total intensity
            sum_out = out.sum()#numpy.sum(numpy.sum(out))
            if sum_out>0:
                out*=flux/sum_out
            
            #scaled_out = out/(float(sum_out)+1E-10)*flux
        else:
            out *= 0#numpy.zeros((pixels,pixels),numpy.float64)
        return out

    def offset_sinc(self,n,d,x,y,fn):
        #This gives the correct fwhm
        #d = fwhm value
        r = d/1.4
        if r<1:
            r = 1.
        #fn = numpy.zeros((n,n),numpy.float64)
        for i in range(n):
            for j in range(n):
                a = (i-x-n/2.+0.5)*numpy.pi
                b = (j-y-n/2.+0.5)*numpy.pi
                if (x==n/2.):
                    c = 1.
                else:
                    c = numpy.sin(a/r)/a
                if (y==n/2.):
                    d = 1.
                else:
                    d = numpy.sin(b/r)/b
                fn[i,j]= c*d
                #if i==0:
                    #print j,a,b,c,d
        return fn
                

    def offset_gaussOrig(self,n,d,x,y,fn):
        r = d/2.34 #This sets d to be the 1/e^2 intensity point
        r2=r**2
        invr2=1/(2.*r2)
        #print n
        #fn = numpy.zeros((n,n),numpy.float64)
        for i in range(n):
            for j in range(n):
                rr = (i-x-n/2.+0.5)**2+(j-y-n/2.+0.5)**2
                b=-rr*invr2#/(2.*(r2))
                if(b>-40):
                    fn[i,j] = numpy.exp(b)
                else:
                    fn[i,j] = 1.0E-12
        return fn

    def offset_gauss(self,n,d,x,y,fn):
        r = d/2.34 #This sets d to be the 1/e^2 intensity point
        r2=r**2
        invr2=1/(2.*r2)
        g=numpy.mgrid[:n,:n]
        rr=(g[0]-x-n/2.+0.5)**2+(g[1]-y-n/2.+0.5)**2
        b=-rr*invr2
        fn[:]=numpy.exp(b)
        #print n
        #fn = numpy.zeros((n,n),numpy.float64)
        return fn


    def lgs_position(self):
        # Calculates off-axis distance of LGS in terms of telescope aperture diameter
        theta = self.LGS_theta
        phi = self.LGS_phi*numpy.pi/180.
        alt = self.LGS_alt
        diam = self.tel_diam
        distance = alt*float(theta/206265.) # small angle approximation
        dist_tel = distance/float(diam) # puts the above in units of telescope diameters
        dist_x = dist_tel*numpy.cos(phi)
        dist_y = dist_tel*numpy.sin(phi)
        return dist_x,dist_y

    def wfs_image(self,off=0,focus=0.0,recalc=0,asInt=0,perSubap=0,verbose=0):
        #print 'Started offset ',off
        # Need an oversized image to allow for spots elongated outside individual subapertures
        if type(self.subapImage)!=type(None) and recalc==0:
            return self.subapImage#already generated...
        image = numpy.zeros((self.npup+(self.subap_pix*2),self.npup+(self.subap_pix*2)),numpy.float64)
        global_off = off # This is the number of pixels to offset each wfs frame by
        cent_x,cent_y = self.lgs_position() # This is diameter of pupils to offset lgs by
        #print 'LGS is',sqrt(cent_x**2+cent_y**2)*self.tel_diam, 'metres off-axis at',self.LGS_alt/1000.,'km'
        #print 'This is an offset angle of',self.LGS_theta,'degrees'
        cent_x = cent_x * self.wfs.shape[0]
        cent_y = cent_y * self.wfs.shape[1]
        # Now centres are in terms of pixels arcoss WFS
        #print cent_x, cent_y, self.wfs.shape
        #raw_input('n:')
        #print off,global_off
        output_subaps = numpy.zeros((self.lenslet.shape[0],self.subap_pix*2,self.subap_pix*2),numpy.float64)
        #temp = 0.0
        # Need to calculate a limit for the minimum
        calc_lens_limit = self.lenslet.shape[0]
        coords =numpy.zeros((2),numpy.float64)
        subap_image = numpy.zeros((self.subap_pix*2,self.subap_pix*2),numpy.float64)
        sliceImage=numpy.zeros((self.subap_pix*2,self.subap_pix*2),numpy.float64)
        for i in range(self.lenslet.shape[0]):
            #print i
            #Find the polar coordinates of  the lenslet subaperture from the centre of the primary
            if verbose==1:
                print "%d/%d"%(i,self.lenslet.shape[0])
            if self.lenslet[i,2]>self.subap_calc_limit:#gets enough light...
                centre = self.wfs.shape[0]/2.
                #cent_x,cent_y = self.lgs_position()
                #coords = centre-self.lenslet[i,:2]
                #coords =numpy.zeros((2),numpy.float64)
                #print cent_x,self.lenslet[i,0],cent_y,self.lenslet[i,1]
                coords[0]=cent_x+(centre-self.lenslet[i,0])
                coords[1]=cent_y+(centre-self.lenslet[i,1])
                offset = numpy.sqrt((coords[0]**2)+(coords[1]**2))
                focus_offset = focus*offset/float(centre)*numpy.sqrt(2)
                angle = numpy.arctan2(coords[0],coords[1])#float(self.direction(coords[1],coords[0]))
                offset_dist = (offset/self.wfs.shape[0])*self.tel_diam
                sinangle=numpy.sin(angle)
                cosangle=numpy.cos(angle)
                #if i<20: print coords[0],coords[1],offset,offset_dist
                #if offset_dist>temp:
                #print i,offset_dist
                #    temp = offset_dist
                #print cent_x,self.lenslet[i,0],cent_y,self.lenslet[i,1],offset,centre,offset_dist
                #raw_input('offsets:')
                # Initialise oversized subaperture image
                subap_image*=0# numpy.zeros((self.subap_pix*2,self.subap_pix*2),numpy.float64)
                for j in range(self.LGS_steps):
                    if verbose==2:
                        print "%d/%d, %d/%d"%(i,self.lenslet.shape[0],j,self.LGS_steps)
                    # Now for each lenslet we determine the offset from lenslet centre
                    # of each range gate step that was defined when the class was initialised
                    alt = self.lgs_return_data[j,0]
                    ang = self.ang_offset(alt,self.LGS_alt,offset_dist)
                    pix_offset = ang/float(self.wfs_scale)
                    x = ((pix_offset+focus_offset) * sinangle)+global_off
                    y = ((pix_offset+focus_offset) * cosangle)+global_off
                    if x**2+y**2 <= 2.25*self.subap_pix*self.subap_pix:
                        # Rayleigh backscatter has sin^2 dependence
                        #flux_scale should be 1 if using sodium lgs.
                        # Very close to 1 for LGS return above 100m
                        if self.laserType=="sodium":
                            flux_scale=1.
                        else:
                            flux_scale = (numpy.sin((numpy.pi/2.)-numpy.arctan(offset_dist/float(alt)))**2)*self.lenslet[i,2]
                        self.subap_im(x,y,self.lgs_return_data[j],self.subap_pix*2,sliceImage)*flux_scale
                        subap_image += sliceImage
            # This is probably redundant
            else: subap_image = numpy.zeros((self.subap_pix*2,self.subap_pix*2),numpy.float64)
                #else:
                    #if i == 1: print j,x,y
            #if self.wfs_noise>=1.:
            #    subap_im =self.add_noise(subap_im)
            # Want to determine centroiding accuracy with offset distance
            #cent_x, cent_y = self.centroid(subap_im)
            #centroids.append((cent_x,cent_y,offset,angle))
            #tot_sig = sum(sum(subap_im))
            #FITS.Write(subap_im, 'subap'+`i`+'.fits')
            # Add the subapertures to an output file ready for inputting to model
            output_subaps[i]=subap_image
            pos_x = self.lenslet[i,0]+self.subap_pix
            pos_y = self.lenslet[i,1]+self.subap_pix
            n = (subap_image.shape[0]/2.)
            #print subap_im.shape, image.shape
            image[int(pos_x-n+0.5):int(pos_x+n+0.5),int(pos_y-n+0.5):int(pos_y+n+0.5)]+=subap_image.astype(numpy.float64)
        # Now slice image back to nominal WFS size
        data = image[self.subap_pix:image.shape[0]-self.subap_pix,self.subap_pix:image.shape[0]-self.subap_pix]
        # Write subap images to file
        #util.FITS.Write(output_subaps,'output_subaps.fits')
        if asInt:
            self.subapImage=(data+0.5).astype(numpy.int32)
        else:
            self.subapImage=data
        self.spotPsf=numpy.zeros((self.nsubx,self.nsubx,self.subap_pix,self.subap_pix),numpy.float32)
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                self.spotPsf[i,j]=data[i*self.subap_pix:(i+1)*self.subap_pix,j*self.subap_pix:(j+1)*self.subap_pix].astype(numpy.float32)
        if perSubap:
            return self.spotPsf
        return self.subapImage

    def sumSubaps(self):
        """Compute the total photon count in each subap.
        """
        #if type(self.sig)==type(None):
        s=numpy.zeros((self.nsubx,self.nsubx),numpy.float32)
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                s[i,j]=numpy.sum(numpy.sum(self.spotPsf[i,j]))
        self.sig=s
        return s
# The following are data analysis tools #################################################################################

    def subaps(self,img):
        nsub = self.nsubx
        pix = int((img.shape[0]/float(nsub))+0.5)
        subs = numpy.zeros((nsub*nsub,pix,pix),numpy.float64)
        fov = self.subap_fov
        for subap in range(self.lenslet.shape[0]):
            min_x = int(self.lenslet[subap,0]-(pix/2.)+0.5)
            min_y = int(self.lenslet[subap,1]-(pix/2.)+0.5)
            max_x = int(self.lenslet[subap,0]+(pix/2.)+0.5)
            max_y = int(self.lenslet[subap,1]+(pix/2.)+0.5)
            #try:
            subs[subap] = img[min_x:max_x,min_y:max_y]
            #except:
                #print subs[subap].shape,img[min_x:max_x,min_y:max_y].shape
        return subs

    def offsets(self,steps,pix):
        #print steps,pix
        # steps is number of steps across a pixel
        # pix is number of pixels to step
        n_frames = int(steps*pix)
        output = numpy.zeros((n_frames,self.npup,self.npup),numpy.float64)
        offset_list = numpy.zeros((n_frames),numpy.float64)
        # First generate the small sub-pixel wfs offsets
        #print 'n_frames',n_frames
        for i in range(steps):
            offs = i/float(steps)
            offset_list[i] = offs
            print 'Starting offset',offs
            output[i] = self.wfs_image(offs)
        # Now we can offset the above by entire pixels
        for i in range(1,pix):
            base = i*steps
            temp_offsets = offset_list[0:steps]+i
            offset_list[base:base+steps] = temp_offsets
            #print 'base',base
            temp = output[0:steps,:-i,:-i]
            #print temp.shape
            output[base:base+steps,i:,i:,]=temp
        return output,offset_list

    def load_analysis(self,infile,filename,resolution,offsets=1):
        fov = self.subap_fov
        pix = self.subap_pix
        scale = fov/float(pix) #arcseconds per pixel
        #pix_offset = int((max_offset/float(fov)*pix/sqrt(2))+1)#always round up
        pix_res = 1/float(resolution) #number of pixels per small step
        wfs_data = util.FITS.Read(infile)[1]
        if offsets == 1:
            calib_offsets = self.calc_offsets(wfs_data[0],scale)
        else:
            calib_offsets = numpy.zeros((self.nsubx,self.nsubx),2)
        output = []
        for offset in range(wfs_data.shape[0]):
            output.append(self.analyse_frame(wfs_data[offset],offset*pix_res,scale,calib_offsets))
        self.savefile(output,filename)
            


    def strehl_analysis(self,filename,repeats):
        fov = self.subap_fov
        pix = self.subap_pix
        scale = fov/float(pix) #arcseconds per pixel
        #pix_offset = int((max_offset/float(fov)*pix/sqrt(2))+1)#always round up
        #pix_res = 1/float(resolution) #number of pixels per small step
        wfs_file = 'WFS_'+filename
        wfs_data = util.FITS.Read(wfs_file)[1]
        strehl_file = 'SR_'+filename
        offset_list = util.FITS.Read(strehl_file)[1][:,0]/float(scale*numpy.sqrt(2))
        calib_offsets = []
        for calibration in range(20):
            offsets = self.calc_offsets(self.add_noise(wfs_data[0]),scale)
            calib_offsets.append(offsets)
        calib_offsets = numpy.array(calib_offsets)
        offsets = numpy.sum(calib_offsets)/20.
        output = []
        for run in range(repeats):
            for offset in range(wfs_data.shape[0]):
                data = self.add_noise(wfs_data[offset])
                if offset ==0:
                    util.FITS.Write(data,'WFS_noisy.fits')
                output_data,cent_x,cent_y,frames = self.analyse_frame((data+0.5).astype(numpy.int32),offset_list[offset],scale,offsets)
                strehl = self.get_strehl(cent_x,cent_y,frames)
                output.append((strehl,offset_list[offset]*scale*numpy.sqrt(2)))
                #output.append(self.analyse_frame(data,offset*pix_res,scale,calib_offsets))
        av_sig = numpy.sum(numpy.sum(wfs_data[0])/float(frames))
        header = ['AVG SIG = '+str(av_sig)[:8]]
        output = numpy.array(output)
        util.FITS.Write(output,'SR_noise_'+filename,header)
    

    def load_analysis_noise(self,infile,filename,resolution,repeats,offsets=1):
        fov = self.subap_fov
        pix = self.subap_pix
        scale = fov/float(pix) #arcseconds per pixel
        #pix_offset = int((max_offset/float(fov)*pix/sqrt(2))+1)#always round up
        pix_res = 1/float(resolution) #number of pixels per small step
        wfs_data = util.FITS.Read(infile)[1]
        if offsets == 1:
            calib_offsets = self.calc_offsets(wfs_data[0],scale)
        else:
            calib_offsets = numpy.zeros((self.nsubx,self.nsubx),2)
        output = []
        for run in range(repeats):
            for offset in range(wfs_data.shape[0]):
                data = self.add_noise(wfs_data[offset])
                output_data,cent_x,cent_y,frames = self.analyse_frame((data+0.5).astype(numpy.int32),offset_list[offset],scale,offsets)
                
                #output.append(self.analyse_frame(data,offset*pix_res,scale,calib_offsets))
        #self.savefile(output,filename)        

    def analysis(self,max_offset,resolution,repeats,filename):
        #max offset is maximum number of arcseconds to offset
        #resolution is number of steps to make across a pixel
        #need to convert these both to arcseconds
        fov = self.subap_fov
        pix = self.subap_pix
        scale = fov/float(pix) #arcseconds per pixel
        pix_offset = int((max_offset/float(fov)*pix/numpy.sqrt(2))+1)#always round up
        pix_res = 1/float(resolution) #number of pixels per small step
        #print scale, pix_offset,pix_res
        wfs_data,offset_list = self.offsets(resolution,pix_offset)
        util.FITS.Write(wfs_data,'WFS_'+filename+'.fits')
        # Need to include a calibration by setting centroid offsets to zero for on-axis spot pattern
        offsets = self.calc_offsets(wfs_data[0],scale)
        output = []
        if self.wfs_noise >=1:
            strehls = numpy.zeros((repeats*wfs_data.shape[0],2),numpy.float64)
            # Need to repeat lots of times to get statistical variation
            for run in range(repeats):
                for offset in range(wfs_data.shape[0]):
                    data = self.add_noise(wfs_data[offset])
                    output_data,cent_x,cent_y,frames = self.analyse_frame((data+0.5).astype(numpy.int32),offset_list[offset],scale,offsets) 
                    output.append(output_data)
                    strehl = self.get_strehl(cent_x,cent_y,frames)
                    index = (run*wfs_data.shape[0])+offset
                    strehls[index,1]=strehl
                    strehls[index,0]=offset_list[offset]*scale*numpy.sqrt(2)
        else:
            # Only need one run if no noise is included
            self.wfs_read = 0.
            strehls = numpy.zeros((wfs_data.shape[0],2),numpy.float64)
            for offset in range(wfs_data.shape[0]):
                
                data = self.add_noise(wfs_data[offset])
                #print offset
                output_data,cent_x,cent_y,frames = self.analyse_frame((data+0.5).astype(numpy.int32),offset_list[offset],scale,offsets) 
                output.append(output_data)
                strehl = self.get_strehl(cent_x,cent_y,frames)
                strehls[offset,1]=strehl
                strehls[offset,0]=offset_list[offset]*scale*numpy.sqrt(2)
                #print strehls[offset,:,],offset
        #self.savefile(output,filename)
        #print strehls.shape,cent_x.shape,cent_y.shape
        #print strehls
        util.FITS.Write(strehls,'SR_'+filename+'.fits')
        return wfs_data

    def calc_offsets(self,frame,scale):
        subs = self.subaps(frame)
        output = []
        for i in range(subs.shape[0]):
            centroids = self.centroid(subs[i])
            x = centroids[0]-(self.subap_pix/2.)-0.5 # actual x-pixel offset from subap centre
            y = centroids[1]-(self.subap_pix/2.)-0.5 # actual y-pixel offset from subap centre
            output.append((x,y)) #x,y offsets in pixels from calibrated position
        return numpy.array(output)

    def analyse_frame(self,frame,offset,scale,calib_offsets):
        subs = self.subaps(frame)
        offs = offset*numpy.sqrt(2)
        output = []
        # Only want to analyse lenslets that are 70% or less vignetted
        #print maximum.reduce(self.lenslet[:,2])
        max_area = maximum.reduce(self.lenslet[:,2])*self.max_vignette
        xs = numpy.zeros((subs.shape[0]),numpy.float64)
        ys = numpy.zeros((subs.shape[0]),numpy.float64)
        frames = 0
        for i in range(subs.shape[0]):
            if self.lenslet[i,2]>max_area:
                frames += 1
                centroids = self.centroid(subs[i])
                x = centroids[0]-(self.subap_pix/2.)-0.5-calib_offsets[i,0]-offset # actual x-pixel offset from calibrated position
                y = centroids[1]-(self.subap_pix/2.)-0.5-calib_offsets[i,1]-offset # actual y-pixel offset from calibrated position
                xs[i]=x
                ys[i]=y
                offsets = (numpy.sqrt((x**2)+(y**2)))*scale #offset in arcseconds from calibrated position
                tot_sig = numpy.sum(numpy.sum(subs[i]))
                masked_sig = numpy.sum(numpy.sum(subs[i][1:-1,1:-1]))
                # gives set offset in arcseconds
                # gives set offset in pixels
                # gives measured centroid deviation from pixel centre
                # gives measured offset in arcseconds
                # want to scale offset to output variances
                offsets = offsets/3600./90./float(self.tel_diam)
                output.append((offs*scale,offset,x,y,offsets,tot_sig,masked_sig))
            else: output.append((offs*scale,0.,0.,0.,0.,0.,0.))
        return output,xs,ys,frames

    def savefile(self,datalist,filename):
        handle = file(filename,'w')
        #print len(datalist),len(datalist[0]),datalist[0]
        #print numpy.array(datalist).shape
        util.FITS.Write(numpy.array(datalist),filename+'.fits')
        for run in datalist:
            for line in run:
                for data in line:
                    handle.write(`data`+'\t')
                handle.write('\n')
        handle.close()

    def get_offsets(self,infile,outfile):
        data_in = file(infile,'r')
        data = data_in.readlines()
        data_b = []
        for line in data:
            data_b.append(string.split(line,'\t'))
        out = numpy.zeros((len(data_b),len(data_b[0])-1),numpy.float64)
        for i in range(out.shape[0]):
            for j in range(out.shape[1]):
                out[i,j]=float(data_b[i][j])
        # Now we have an numpy.array of the data
        actual = out[:,0]
        measured = out[:,4]
        # Slice required data
        output = numpy.zeros((actual.shape[0],2),numpy.float64)
        output[:,0]=actual
        output[:,1]=measured

    def read_data(self,infile,outfile):
        data_in = util.FITS.Read(infile)[1]
        actual = data_in[:,:,0] # actual offsets in arcseconds
        measured = data_in[:,:,4] # measured offsets in arcseconds
        out = numpy.zeros((data_in.shape[0],data_in.shape[1]+1),numpy.float64)
        for i in range(data_in.shape[0]):
            out[i,0] = actual[i,0]
            out[i,1:,]= measured[i]
        # Now format for spreadsheetd input
        fp = file(outfile,'w')
        for line in out:
            for el in line:
                fp.write(`el`+'\t')
            fp.write('\n')
        fp.close()

    def get_strehl(self,cent_x,cent_y,frames):
        pix_scale = self.subap_fov/self.subap_pix #arcseconds per pixel
        sorpist=numpy.zeros((30,30),numpy.float64)
        offx=reshape(cent_x,(30,30))*pix_scale*0.59#-cent_x[0]-pix_actual[i]*sqrt(2.),(32,32))
        offy=reshape(cent_y,(30,30))*pix_scale*0.59#-cent_y[0]-pix_actual[i]*sqrt(2.),(32,32))

        sor30.fit(offx,offy,sorpist)

        phsvar = numpy.sum(numpy.sum(sorpist*sorpist)) / float(frames)

        if phsvar==0:
            strehl = numpy.exp(-phsvar)
        else:
            strehl = numpy.sqrt(float(1+(phsvar**2)))-phsvar

        return strehl


    def plot_cents(self,filename):
        #print filename
        data = util.FITS.Read(filename)[1]
        actual = data[:,:,0] #arcseconds
        offs = data[:,:,4] #arcseconds
        pix_actual = data[:,:,1] #pixels
        cent_x = data[:,:,2] #pixels
        cent_y = data[:,:,3] # pixels

        pix_scale = self.subap_fov/self.subap_pix #arcseconds per pixel
                
        for i in range(0,len(cent_x)):
            sorpist=numpy.zeros((30,30),numpy.float64)
            #offx=reshape(offs[i],(32,32)) / sqrt(2.)
            #offy=reshape(offs[i],(32,32)) / sqrt(2.)
            offx=reshape(cent_x[i],(30,30))*pix_scale*0.59#-cent_x[0]-pix_actual[i]*sqrt(2.),(32,32))
            offy=reshape(cent_y[i],(30,30))*pix_scale*0.59#-cent_y[0]-pix_actual[i]*sqrt(2.),(32,32))
            #print offx

            sor30.fit(offx,offy,sorpist)
            
            window(0,wait=1)
            pli(sorpist)
            
            phsvar = numpy.sum(numpy.sum(sorpist*sorpist)) / 684.
            try:
                print i,numpy.exp(-phsvar)
            except:
                print 0.0

            
            #raw_input('N:')


class Flux:
    """from photon_return.py - TJM
    Calculates the flux from LGS.
    """
    def __init__(self,wavelength=514.,tel_trans=0.67,optics_trans=0.44,power=30.,launch_trans=0.8,pulse_rate=5000.,frame_rate=300.,QE=0.87,zenith=0.,tel_alt=2269.,scale_height=7400.,sodProfileResolution=10,sodProfilePeak1=90000.,sodProfilePeak2=100000.,sodProfileWidth1=3000.,sodProfileWidth2=3000.,sodProfileStrength1=0.4,sodProfileStrength2=0.6,numericProfile=None):
        """
        tel_alt is the altitude of the telescope (ie usually on top of a mountain somewhere!)
        wavelength is the wavelength to use, nm
        tel_trans is telescope transmission function
        optics_trans is optics transmission function
        power is laser power in W.
        pulse_rate is laser pulse rate in Hz
        launch_trans is launch optics transmission function
        frame_rate is WFS frame rate.
        QE is WFS CCD QE (0 to 1)
        zenith is the pointing angle of telescope (degrees)
        scale_height is atmospheric scale height (decay of pressure with altitude) (ask TJM)
        numericProfile is optional, a LGS sodium profile.  This is 2 arrays, of (heights, relative strengths).


        For sodium, the following parameters affect result:
        tel_trans, optics_trans, power, launch_trans, QE, frame_rate, zenith, tel_alt, *Profile*.
        """
        self.wavelength=wavelength
        self.tel_trans=tel_trans
        self.optics_trans=optics_trans
        self.power=power
        self.pulse_rate=pulse_rate
        self.launch_trans=launch_trans
        self.pulse_rate=pulse_rate
        self.frame_rate=frame_rate
        self.QE=QE
        self.zenith=zenith
        self.tel_alt=tel_alt
        self.scale_height=scale_height
        self.sodProfileResolution=sodProfileResolution
        self.sodProfilePeak1=sodProfilePeak1
        self.sodProfilePeak2=sodProfilePeak2
        self.sodProfileWidth1=sodProfileWidth1
        self.sodProfileWidth2=sodProfileWidth2
        self.sodProfileStrength1=sodProfileStrength1
        self.sodProfileStrength2=sodProfileStrength2
        self.numericProfile=numericProfile
        self.calcSodiumFlux()
        h = 6.626068E-34 # Planck's Constant
        c = 3E8     # Speed of light

        L = float(self.wavelength * 1.E-09)
        # Look at Bucholtz Applied Optics v34 n15 1995
        if L <=500E-09:
            A = 3.01577E-028
            B = 3.55212
            C = 1.35579
            D = 0.11563
        else:
            A = 4.01061E-028
            B = 3.99668
            C = 1.10298E-03
            D = 2.71393E-02
        Lm = L * 1E6 # wavelength in microns
        scatter_csa = A*(Lm**(-B+(C*Lm)+(D/Lm)))/(100.**2) #converts to m^2
        #print scatter_csa, 5E-57*(L**-4)

        tel_trans = self.tel_trans
        opt_trans = self.optics_trans

        per_pulse = (self.power/float(self.pulse_rate))*self.launch_trans

        self.pulses_per_frame = self.pulse_rate/float(self.frame_rate)

        #Half the LIDAR equation is independent of altitude
        self.LIDAR_base = (per_pulse*scatter_csa*tel_trans*opt_trans*self.QE*L)/(h*c*4.*numpy.pi)

        # This is used to scale the scale height of the atmosphere
        # The scale height of the atmosphere is the height at which
        # density reaches 1/e^2 of its sea-level value
        self.zen = 1/numpy.cos(self.zenith*numpy.pi/float(180.))

        sea_level_dens = 1.2928 #kg/m^3
        Na = 6.023E23 # Avogadro's number
        molecular_mass = 28.97E-03 #kg/mol

        self.base_dens = sea_level_dens*Na/float(molecular_mass)

        # The following is taken from La Palma technical note 31
        # Referencing Hayes and Latham (1975) Ap J 197, 593
        l = (self.wavelength/1000.)**-2.
        n = (0.23465+(1.076E2/(146-l))+(0.93161/(41-l)))**2
        self.base_trans = 9.4977E-3*(l**2)*n

        self.T_oz = 0.012550 # Mag per airmass extinction due to ozone at 510nm
        # Taken from Gast (1960) Handbook of Geophysics USAF
        # and Allen (1963) Astrophysical Quantities
        T_obs = (self.base_trans*numpy.exp(-self.tel_alt/float(self.scale_height)))+self.T_oz
        self.base_extinction = 100.**(-T_obs/5.)

    def gamma(self):
        # This is a function of wavelength
        # Look at Bucholtz Applied Optics v34 n15 1995
        # We're only using green so
        return 1.442E-02

    def phase_dispersal(self,angle):
        A = 3/4.*float(1+(2*self.gamma))
        B = 1 + (3*self.gamma)
        C = (1-self.gamma)
        D = numpy.cos(angle)**2
        return A*(B+(C*D))

    def mol_dens(self,base_alt,alt):
        """Molar density at a given height?"""
        a = base_alt+alt
        out = self.base_dens/numpy.exp(a/float(self.scale_height*self.zen))
        return out #,2.6E29/(a*exp(a/self.scale_height))

    def atmos_trans(self,base_alt,alt):
        """Atmospheric transmission?"""
        T_alt = (self.base_trans*numpy.exp(-(alt+base_alt)/float(self.scale_height)))+self.T_oz
        extinction = 100**(-T_alt/5.)
        transmission = 1-(extinction-self.base_extinction)
        return transmission
        #if alt < 10000.:
        #    return 0.9
        #else:
        #    return 0.8 #base transmission for 20km LGS

    def rayleigh(self,slice_alt,slice_depth):
        """Return photons per m^2 per frame"""
        dens = self.mol_dens(self.tel_alt,slice_alt)
        trans = self.atmos_trans(self.tel_alt,slice_alt)**2
        return_per_pulse = self.LIDAR_base*dens*trans*slice_depth/float(slice_alt**2)
        # Need to multiply this value by the area of the subaperture in m^2
        # Also need the angular offset as backscatter has an angular component
        # equal to sin^2(off-axis angle)
        # These values are calculated in the LGS WFS program
        out = return_per_pulse*self.pulses_per_frame
        print "Rayleigh %g %g %g"%(slice_alt,slice_depth,out)
        return out

    def calcSodiumFlux(self):
        """
        Taken from book Optics in Astrophysics, Renaud Foy, page 220 (Laser tech for LG AO chapter):
        Equation is Flux (photons/m^2/s) = n c sig Pl cos(z) / ( 4pi H^2 E)
        where n is the overall efficiency factor (0.1), c is 5e13/m^2, sig is the peak D2 cross section (8.8e-16), Pl is the laser power (Watts), z is the zenith angle, H is the LGS altitude in m, and E is the photon enerhgy, 3.38e-19 J.
        Note, these values were taken from the book.  So, using this, we can get the total flux.  Then, we divide this flux over the number of slices in the sodium profile, which then gives us the flux at this slice alt, for this slice depth.
        """
        n=self.tel_trans*self.optics_trans*self.launch_trans*self.QE#efficiency factor
        C=5e13#m^-2
        sig=8.8e-16#m^2
        power=self.power#watts
        zenithAngle=self.zenith#degrees
        Ephot=3.38e-19#joules
        height=92000.-self.tel_alt#m. - height of the sodium layer...
        self.totSodiumFlux=n*C*sig*power*numpy.cos(zenithAngle*numpy.pi/180.)/(4*numpy.pi*height**2*Ephot)        
    def sodium(self,sliceAltList):
        """Return the number of photons per m^2 per frame.
        First, get the sodium density at listed heights - 2 gausians with peaks at 85 and 95km.
        Then work out atmospheric transmission.
        Work out fraction of excited from this.
        Then scale with 1/r^2.
        Taken from book Optics in Astrophysics, Renaud Foy, page 220 (Laser tech for LG AO chapter):
        Equation is Flux (photons/m^2/s) = n c sig Pl cos(z) / ( 4pi H^2 E)
        where n is the overall efficiency factor (0.1), c is 5e13/m^2, sig is the peak D2 cross section (8.8e-16), Pl is the laser power (Watts), z is the zenith angle, H is the LGS altitude in m, and E is the photon enerhgy, 3.38e-19 J.
        Note, these values were taken from the book.  So, using this, we can get the total flux.  Then, we divide this flux over the number of slices in the sodium profile, which then gives us the flux at this slice alt, for this slice depth.
        
        """
        #trans = self.atmos_trans(self.tel_alt)**2

        #Now find out how much of this flux is attributed to each slice.
        profile=[]
        #pad the slices - so that we can get slice width of the start and ending slices...
        sliceAltList=[sliceAltList[0]*2-sliceAltList[1]]+sliceAltList+[sliceAltList[-1]*2-sliceAltList[-2]]
        for i in range(1,len(sliceAltList)-1):
            sliceAlt=sliceAltList[i]
            sliceDepth=(sliceAltList[i+1]-sliceAltList[i-1])/2.
            fr=(sliceAlt+sliceAltList[i-1])/2.
            to=(sliceAlt+sliceAltList[i+1])/2.
            profile.append(self.integrateSodiumProfile(fr,to))
        profile=numpy.array(profile).astype(numpy.float64)
        profile/=numpy.sum(profile)
        profile*=self.totSodiumFlux
        return profile/self.frame_rate
    def integrateSodiumProfile(self,fr,to):
        """Here, we integrated the sodium profile (double gaussian) from fr to to.
        """
        resolution=self.sodProfileResolution#number of numerical integration points to use...
        rnge=to-fr
        step=rnge/resolution
        x=fr+step/2.
        if self.numericProfile==None:
            integ=0.
            for i in range(resolution):
                integ+=self.sodProfileStrength1*numpy.exp(-((x-self.sodProfilePeak1)/self.sodProfileWidth1)**2/2)
                integ+=self.sodProfileStrength2*numpy.exp(-((x-self.sodProfilePeak2)/self.sodProfileWidth2)**2/2)
                x+=step
        else:#Use the numeric profile.
            
            x=self.numericProfile[0]
            y=self.numericProfile[1]
            interp=scipy.interpolate.interp1d(x,y,"cubic",bounds_error=False,fill_value=0.)
            xx=numpy.arange(resolution)*step+fr
            yy=interp(xx)
            yy=numpy.where(yy<0,0,yy)
            integ=yy.sum()
            
        return integ


    def plotSodiumProfile(self,start=0,height=120000):
        height=int(height)
        start=int(start)
        a=numpy.zeros((height-start,),"d")
        if self.numericProfile==None:
            for x in range(start,height):
                #print x
                a[x-start] =self.sodProfileStrength1*numpy.exp(-((x-self.sodProfilePeak1)/self.sodProfileWidth1)**2/2)
                a[x-start]+=self.sodProfileStrength2*numpy.exp(-((x-self.sodProfilePeak2)/self.sodProfileWidth2)**2/2)
        else:
            x=self.numericProfile[0]
            y=self.numericProfile[1]
            interp=scipy.interpolate.interp1d(x,y,"cubic",bounds_error=False,fill_value=0.)
            xx=range(start,height)
            a=interp(xx)
            a=numpy.where(a<0,0,a)
        return a
        
if __name__=="__main__":
    import sys
    gsList=None
    layerList=None
    telDiam=None
    telSec=None
    fill=False
    if len(sys.argv)>1:
        gsList=eval(sys.argv[1])
    if len(sys.argv)>2:
        layerList=eval(sys.argv[2])
    if len(sys.argv)>3:
        telDiam=eval(sys.argv[3])
    if len(sys.argv)>4:
        telSec=eval(sys.argv[4])
    if len(sys.argv)>1:
        fill=eval(sys.argv[5])
    displayGSOverlap(gsList,layerList,telDiam,telSec,fill)
