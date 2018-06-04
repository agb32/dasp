import base.readConfig
base.readConfig.init(globals())
#For physical optics propagation, the phase screen needs to be larger than the pupil, and needs to be sampled at approx 1cm or so.
this.infScrn=new()
this.infScrn.npup=400
this.infScrn.telDiam=4.2
this.physProp=new()
this.physProp.npup=400
this.physProp.npupClipped=100
this.physProp.telDiam=4.2
this.physProp.ntel=400

#The rest of the simulation sees a pupil of 100x100 phase elements for a 1.05m telescope.
npup=100
telDiam=4.2/4.
telSec=telDiam/7.#Central obscuration
ntel=npup#Telescope diameter in pixels
ngsLam=640.#NGS wavelength
sciLam=640.#sci wavelength - for fresnel, probably has to equal ngsLam unless a different atmos module is used.
import util.tel
pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam)

tstep=1/250.#Simulation timestep in seconds (250Hz).
AOExpTime=40.#40 seconds exposure (use --iterations=xxx to modify)


from util.atmos import geom,layer,source
atmosDict={}
nlayer=2 #2 atmospheric layers
layerList={"allLayers":["L%d"%x for x in range(nlayer)]}
strList=[0.5,0.5]#+[0.1]*(nlayer-1)#relative strength of the layers
hList=[4000.,10000.]#range(0,nlayer*1000,1000)#height of the layers
vList=[10.,13.]
dirList=[0.,30.]
for i in range(nlayer):
 atmosDict["L%d"%i]=layer(hList[i],dirList[i],vList[i],strList[i],10+i)

l0=10. #outer scale
r0=0.137 #fried's parameter

atmosPhaseType="phaseamp"


#Create the WFS overview
wfs_nsubx=5 #Number of subaps
import util.guideStar
wfsDict={"1":util.guideStar.NGS("1",wfs_nsubx,0.,0.,phasesize=npup/wfs_nsubx,nimg=8,minarea=0.5,sig=2e4,sourcelam=ngsLam,reconList=["recon"],pupil=pupil,atmosPhaseType="phaseamp")}
wfsOverview=util.guideStar.wfsOverview(wfsDict)

#Create a Science overview.
import util.sci
sciDict={}
nsci=this.getVal("nsci",1)
if nsci==1:
 phslam=ngsLam
else:
 phslam=sciLam
for i in range(nsci):
 sciDict["sci%d"%(i+1)]=util.sci.sciInfo("sci%d"%(i+1),i*10.,0.,pupil,sciLam,phslam=phslam,phaseType="phaseamp")
 sciOverview=util.sci.sciOverview(sciDict)


sourceList=[]
#the wfs
sourceList.append(wfsOverview.getWfsByID("1"))

#and psf
for i in range(nsci):
 sourceList.append(sciOverview.getSciByID("sci%d"%(i+1)))


atmosGeom=geom(atmosDict,sourceList,ntel,npup,telDiam,r0,l0)
this.physProp.atmosGeom=geom(atmosDict,sourceList,this.physProp.ntel,this.physProp.npup,this.physProp.telDiam,r0,l0,ignoreZenithWarning=1)
this.infScrn.atmosGeom=this.physProp.atmosGeom



#Create the DM object.
nAct=wfs_nsubx+1
from util.dm import dmOverview,dmInfo
dmInfoList=[dmInfo('dm',[x.idstr for x in sourceList],0.,nAct,minarea=0.1,actuatorsFrom="recon",pokeSpacing=(None if wfs_nsubx<20 else 10),maxActDist=1.5,decayFactor=0.95)]
dmOverview=dmOverview(dmInfoList,atmosGeom)

rcond=0.05#condtioning value for SVD
recontype="pinv"#reconstruction type
pokeval=1.#strength of poke
gainFactor=0.5#Loop gain
computeControl=1#To compute the control matrix after poking
reconmxFilename="rmx.fits"#control matrix name (will be created)
pmxFilename="pmx.fits"#interation matrix name (will be created)

