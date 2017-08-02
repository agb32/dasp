import base.readConfig
base.readConfig.init(globals())
#For physical optics propagation, the phase screen needs to be larger than the pupil, and needs to be sampled at approx 1cm or so.
this.infScrn=new()
this.infScrn.npup=160
this.infScrn.telDiam=0.8
this.physProp=new()
this.physProp.npup=160
this.physProp.npupClipped=40
this.physProp.telDiam=0.8
this.physProp.ntel=160

#The rest of the simulation sees a pupil of 100x100 phase elements for a 1.05m telescope.
npup=40
telDiam=0.8/4.
telSec=telDiam/7.#Central obscuration
ntel=npup#Telescope diameter in pixels
sciLam=640.#sci wavelength
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


#Create a Science overview.
import util.sci
sciDict={}
nsci=this.getVal("nsci",1)

for i in range(nsci):
 sciposx=i*10.
 sciposy=0.
 sciposr=numpy.sqrt(sciposx**2+sciposy**2)
 sciposphi=numpy.arctan2(sciposy,sciposx)*180/numpy.pi
 sciDict["sci%d"%(i+1)]=util.sci.sciInfo("sci%d"%(i+1),sciposr,sciposphi,pupil,sciLam,phslam=sciLam,phaseType="phaseamp",inboxDiamList=[],psfFilename="psf%d.fits"%(i+1))
 sciOverview=util.sci.sciOverview(sciDict)


sourceList=[]

#and psf
for i in range(nsci):
 sourceList.append(sciOverview.getSciByID("sci%d"%(i+1)))


atmosGeom=geom(atmosDict,sourceList,ntel,npup,telDiam,r0,l0)
this.physProp.atmosGeom=geom(atmosDict,sourceList,this.physProp.ntel,this.physProp.npup,this.physProp.telDiam,r0,l0,ignoreZenithWarning=1)
this.infScrn.atmosGeom=this.physProp.atmosGeom



