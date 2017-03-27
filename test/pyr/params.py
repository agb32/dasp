
import base.readConfig
base.readConfig.init(globals())
wfs_nsubx=10 #Number of subaps
tstep=1/250.#Simulation timestep in seconds (250Hz).
AOExpTime=40.#40 seconds exposure (use --iterations=xxx to modify)
npup=wfs_nsubx*8#Number of phase points across the pupil
telDiam=4.2
telSec=telDiam/7.#Central obscuration
ntel=npup#Telescope diameter in pixels
nAct=wfs_nsubx+1#Number of actuators across the DM
ngsLam=640.#NGS wavelength
sciLam=1650.#Science wavelength
nsci=1
import util.tel
#Create a pupil function
pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam)

if hasattr(this.globals,"pyrSteps"):
 pyrSteps=this.globals.pyrSteps
else:
 pyrSteps=8
if hasattr(this.globals,"pyrModAmp"):
 pyrModAmp=this.globals.pyrModAmp
else:
 pyrModAmp=8.

#Create the WFS overview
import util.guideStar
wfsDict={"1":util.guideStar.NGS("1",wfs_nsubx,0.,0.,phasesize=npup/wfs_nsubx,nfft=npup*2,clipsize=npup*2,nimg=wfs_nsubx*2,                                minarea=0.5,sig=1e6,sourcelam=ngsLam,                                reconList=["recon"],pupil=pupil,pyrSteps=pyrSteps,pyrModAmp=pyrModAmp)}
wfsOverview=util.guideStar.wfsOverview(wfsDict)

#Create a Science overview.
import util.sci
sciDict={}
if nsci==1:
 phslam=ngsLam
else:
 phslam=sciLam
for i in range(nsci):
 sciDict["sci%d"%(i+1)]=util.sci.sciInfo("sci%d"%(i+1),i*10.,0.,pupil,sciLam,phslam=phslam)
 sciOverview=util.sci.sciOverview(sciDict)

#Create the atmosphere object and source directions.
from util.atmos import geom,layer,source
atmosDict={}
nlayer=2 #2 atmospheric layers
layerList={"allLayers":["L%d"%x for x in range(nlayer)]}
strList=[0.9]+[0.1]*(nlayer-1)#relative strength of the layers
hList=range(0,nlayer*1000,1000)#height of the layers
vList=[10.]*nlayer#velocity of the layers
dirList=[0.]*nlayer#direction (degrees) of the layers
for i in range(nlayer):
 atmosDict["L%d"%i]=layer(hList[i],dirList[i],vList[i],strList[i],10+i)
sourceList=[]
#the wfs
sourceList.append(wfsOverview.getWfsByID("1"))

#and psf
for i in range(nsci):
 sourceList.append(sciOverview.getSciByID("sci%d"%(i+1)))
l0=10. #outer scale
r0=0.137 #fried's parameter
atmosGeom=geom(atmosDict,sourceList,ntel,npup,telDiam,r0,l0)


#Create the DM object.
from util.dm import dmOverview,dmInfo
dmInfoList=[dmInfo('dm',[x.idstr for x in sourceList],0.,nAct,minarea=0.1,actuatorsFrom="recon",                   pokeSpacing=(None if wfs_nsubx<20 else 10),maxActDist=1.5,decayFactor=0.95)]
dmOverview=dmOverview(dmInfoList,atmosGeom)

#reconstructor
this.tomoRecon=new()
r=this.tomoRecon
r.rcond=0.05#condtioning value for SVD
r.recontype="pinv"#reconstruction type
r.pokeval=1.#strength of poke
if hasattr(this.globals,"gain"):
 r.gainFactor=this.globals.gain
else:
 r.gainFactor=2.0#Loop gain - of 2! Seems better than 0.5, probably conditioning
r.computeControl=1#To compute the control matrix after poking
r.reconmxFilename="rmx.fits"#control matrix name (will be created)
r.pmxFilename="pmx.fits"#interation matrix name (will be created)
