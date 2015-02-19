useCmod=1
wfs_nsubx=6
tstep=1/250.
l0=10.
r0=0.137
wfs_sig=1e6
npup=wfs_nsubx*8
telDiam=4.2*wfs_nsubx/6.
ntel=npup
wfs_n=npup/wfs_nsubx
wfs_nfft=wfs_n*2
wfs_clipsize=wfs_nfft
wfs_nimg=wfs_clipsize/2
wfs_ncen=wfs_nimg
nAct=wfs_nsubx+1
ngsLam=640.
nlayer=1
ndm=nlayer
strList=[1.]
hList=[0.]
vList=[10.]
dirList=[0.]
sciLam=1650.
if wfs_nsubx<20:
 pokeSpacing=None
else:
 pokeSpacing=10
from util.atmos import geom,layer,source
d={}
for i in range(nlayer):
 d["L%d"%i]=layer(hList[i],dirList[i],vList[i],strList[i],10+i)
sourceList=[]
#the wfs
sourceList.append(source("a",0.,0.,-1,wfs_nsubx,ngsLam,reconList=["ngs"]))
#and psf
sourceList.append(source("m",0.,0.,-1,None,sciLam,phslam=ngsLam))
atmosGeom=geom(d,sourceList,ntel,npup,telDiam,r0,l0)



from util.dm import dmOverview,dmInfo,calcActuators
all=atmosGeom.sourceDict.keys()
all.sort()
reconList=["a"]
dmInfoList=[]
dmInfoList.append(dmInfo('dm',all,0.,nAct,minarea=0.1,closedLoop=1,actSpacing=None,reconLam=ngsLam,actuatorsFrom=["ngs"],\
	reconstructList=reconList,pokeSpacing=pokeSpacing,interpType="spline",maxActDist=1.5,slaving="auto"))
dmObj=dmOverview(dmInfoList,atmosGeom)


telSec=telDiam/7.
scinSamp=10
AOExpTime=40.
simulationTime=0.
wfs_minarea=0.5
import util.tel
pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam,spider=None)

scrnScale=telDiam/float(ntel)


calsource=0
timing=0
wfs_int=tstep
wfs_lat=0.
wfs_read_mean=0.
wfs_read_sigma=0
wfs_floor=wfs_read_mean+0*wfs_read_sigma
wfs_skybrightness=0.
nthreads="all"
this.infAtmos=new()
this.infAtmos.dataType="d"
this.infScrn=new()
this.infScrn.dataType="d"
this.infScrn.seed=None
this.iscrn=new()
this.iscrn.seed=None

this.science=new()
s=this.science
s.simFilename=None
s.science_integrate=1
s.science_calcRMS=0
s.zero_science=10
s.hist_list_size=int(AOExpTime/tstep/scinSamp)
s.fitsFilename=None
s.scicsvFilename="results.csv"
s.scinfft=npup*2
s.scinimg=s.scinfft
s.inboxDiamList=[0.1]
s.usedmpup=0
s.sourceLam=sciLam
s.histFilename="resultshist.fits"

this.tomoRecon=new()
r=this.tomoRecon
r.rcond=0.05
r.recontype="pinv"
r.pokeval=1.
r.gainFactor=0.5
r.decayFactorOpen=0.
r.computeControl=1
r.dmModeType="poke"
r.reconmxFilename="pmx.fits"
r.pmxFilename="rms.fits"
r.abortAfterPoke=0
