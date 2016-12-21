wfs_nsubx=6 #Number of subaps
tstep=1/250.#Simulation timestep in seconds (250Hz).
AOExpTime=40.#40 seconds exposure (use --iterations=xxx to modify)
npup=wfs_nsubx*8#Number of phase points across the pupil
telDiam=4.2*wfs_nsubx/6.#Telescope diameter
telSec=telDiam/7.#Central obscuration
ntel=npup#Telescope diameter in pixels
nAct=wfs_nsubx+1#Number of actuators across the DM
ngsLam=640.#NGS wavelength
sciLam=1650.#Science wavelength
import util.tel
#Create a pupil function
pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam)

#Create the WFS overview
import util.guideStar
wfsDict={"a":util.guideStar.NGS("a",wfs_nsubx,0.,0.,phasesize=npup/wfs_nsubx,\
                                minarea=0.5,sig=1e6,sourcelam=ngsLam,\
                                reconList=["ngs"],pupil=pupil)}
wfsOverview=util.guideStar.wfsOverview(wfsDict)

#Create a Science overview.
import util.sci
sciDict={"m":util.sci.sciInfo("m",0.,0.,pupil,sciLam,phslam=ngsLam)}
sciOverview=util.sci.sciOverview(sciDict)

#Create the atmosphere object and source directions.
from util.atmos import geom,layer,source
atmosDict={}
nlayer=1 #1 atmospheric layer
strList=[1.]#relative strength of the layers
hList=[0.]#height of the layers
vList=[10.]#velocity of the layers
dirList=[0.]#direction (degrees) of the layers
for i in range(nlayer):
 atmosDict["L%d"%i]=layer(hList[i],dirList[i],vList[i],strList[i],10+i)
sourceList=[]
#the wfs
sourceList.append(wfsOverview.getWfsByID("a"))
#and psf
sourceList.append(sciOverview.getSciByID("m"))
l0=10. #outer scale
r0=0.137 #fried's parameter
atmosGeom=geom(atmosDict,sourceList,ntel,npup,telDiam,r0,l0)


#Create the DM object.
from util.dm import dmOverview,dmInfo
dmInfoList=[dmInfo('dm',['a','m'],0.,nAct,minarea=0.1,actuatorsFrom="ngs",\
                   pokeSpacing=(None if wfs_nsubx<20 else 10),maxActDist=1.5,interpType="hex")]
dmObj=dmOverview(dmInfoList,atmosGeom)

seed=1

#reconstructor
this.tomoRecon=new()
r=this.tomoRecon
r.rcond=0.05#condtioning value for SVD
r.recontype="pinv"#reconstruction type
r.pokeval=1.#strength of poke
r.gainFactor=0.5#Loop gain
r.computeControl=1#To compute the control matrix after poking
r.reconmxFilename="rmx.fits"#control matrix name (will be created)
r.pmxFilename="pmx.fits"#interation matrix name (will be created)
