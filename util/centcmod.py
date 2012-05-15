import cmod.cent
import numpy

CALSOURCE=1
SIG=2
ADDPOISSON=3
NCEN=4
READNOISE=5
READBG=6
NOISEFLOOR=7
SKYBRIGHTNESS=8
PXLPOWER=9
SPOTPSF=10
NTHREADS=11
OPTICALBINNING=12
CENTWEIGHT=13
CORRELATIONCENTROIDING=14
CORRTHRESH=15
CORRPATTERN=16
CALDATA=17
REFCENTS=18
CALCOEFF=19
USEBRIGHTEST=20
class centcmod:
  def __init__(self,nthreads,nsubx,ncen,fftsize,clipsize,
               nimg,phasesize,readnoise,readbg,
               addPoisson,noiseFloor,sig,skybrightness,
               calsource,pxlPower,nintegrations,seed,
               phs,pup,spotpsf,cents,bimg,minarea,opticalBinning,centWeight,correlationCentroiding,corrThresh,corrPattern,corrimg,threshType,imageOnly,useBrightest):
    """Wrapper for the c centroid module.
    Here, sig can be a float or a float array.
    """
    self.nthreads=nthreads
    self.nsubx=nsubx
    self.nsubaps=nsubx*nsubx
    self.ncen=ncen
    self.fftsize=fftsize
    self.clipsize=clipsize
    self.nimg=nimg
    self.phasesize=phasesize
    self.readnoise=readnoise
    self.readbg=readbg
    self.addPoisson=addPoisson
    self.noiseFloor=noiseFloor
    self.sig=sig
    self.skybrightness=skybrightness
    self.calsource=calsource
    self.minarea=minarea
    self.opticalBinning=opticalBinning
    if pxlPower==None:
      pxlPower=1.
    self.pxlPower=pxlPower
    
    self.nintegrations=nintegrations
    self.seed=seed
    self.phs=phs
    
    self.pup=pup
    self.pupfn=pup.perSubap(nsubx,vectorAlign=0)
    self.subflag=pup.getSubapFlag(nsubx,minarea).astype(numpy.int32)
    self.fracSubArea=(pup.getSubarea(nsubx)/(phasesize*phasesize)).astype(numpy.float32)
    self.spotpsf=spotpsf
    self.cents=cents
    self.bimg=bimg
    if type(centWeight)==type(""):#will be created later...
      centWeight=None
    self.centWeight=centWeight
    self.correlationCentroiding=correlationCentroiding
    self.corrThresh=corrThresh
    self.corrPattern=corrPattern
    self.corrimg=corrimg
    self.threshType=threshType
    self.imageOnly=imageOnly
    self.useBrightest=useBrightest
    if correlationCentroiding:
      if corrPattern!=None:
        if corrPattern.dtype.char!="f" or corrPattern.shape!=(nsubx,nsubx,nimg,nimg):
          raise Exception("corrPattern wrong shape or type")
      if corrimg.dtype.char!="f" or corrimg.shape!=(nsubx,nsubx,nimg,nimg):
        raise Exception("corrimg wrong shape or type")

    print "centcmod: cmod.cent.initialise"
    self.centstruct=cmod.cent.initialise(self.nthreads,self.nsubaps,self.ncen,self.fftsize,self.clipsize,self.nimg,self.phasesize,self.readnoise,self.readbg,self.addPoisson,self.noiseFloor,self.sig,self.skybrightness,self.calsource,self.pxlPower,self.nintegrations,self.seed,self.phs,self.pupfn,self.spotpsf,self.cents,self.subflag,self.bimg,self.fracSubArea,self.opticalBinning,self.centWeight,correlationCentroiding,corrThresh,corrPattern,corrimg,threshType,imageOnly,useBrightest)

  def run(self,calsource):
    if calsource!=self.calsource:
      cmod.cent.update(self.centstruct,CALSOURCE,calsource)
      self.calsource=calsource
    t=cmod.cent.run(self.centstruct)
    return t
  def free(self):
    if self.centstruct!=None:
      cmod.cent.free(self.centstruct)
      self.centstruct=None
  def update(self,what,val):
    """what should be an int as defined at the top of this value, eg CALSOURCE=1, SIG=2, etc.  
    val is the new value.
    Note, this is not for general use - does not check that what is legal, and does not update the python values, meaning cmod and python becomes out of sync.  
    Currently, only used when doing SHS calibration.
    """
    cmod.cent.update(self.centstruct,what,val)

if __name__=="__main__":
  import Numeric
  import util.centcmod
  import util.tel
  import gist
  nthreads=2
  nsubx=8
  nsubaps=nsubx**2
  ncen=6
  fftsize=32
  nimg=32
  phasesize=8
  readnoise=1.
  readbg=1.
  addPoisson=1
  noiseFloor=1.
  sig=1000.
  skybrightness=1.
  calsource=1
  pxlPower=1.
  nintegrations=1
  seed=1
  pup=util.tel.Pupil(nsubx*phasesize,nsubx*phasesize/2,0,nsubx)
  #subflag=pup.subflag.astype(Numeric.Int32)
  #fracSubArea=(pup.subarea/phasesize**2).astype(Numeric.Float32)
  #pupfn=pup.perSubap(nsubx)
  cents=Numeric.zeros((nsubx*nsubx,2),Numeric.Float32)
  phs=Numeric.zeros((nsubx*phasesize,nsubx*phasesize),Numeric.Float32)
  phs[:,]=(Numeric.arange(nsubx*phasesize)*10./15.).astype("f")
  reorderedPhs=Numeric.zeros((nsubx,nsubx,nintegrations,phasesize,(phasesize+3)&~3),Numeric.Float32)
  for i in range(nsubx):
    for j in range(nsubx):
      for k in range(nintegrations):
        reorderedPhs[i,j,k,:,:phasesize]=phs[i*phasesize:(i+1)*phasesize,j*phasesize:(j+1)*phasesize]


  spotpsfdim=4
  if spotpsfdim==2:
    spotpsf=Numeric.zeros((fftsize,fftsize),"f")
    spotpsf[fftsize/2-5:fftsize/2+5,fftsize/2-5:fftsize/2+5]=1
    spotpsf[fftsize/2-7:fftsize/2-3,fftsize/2-7:fftsize/2-3]=1
    spotpsf[fftsize/2+3:fftsize/2+7,fftsize/2+3:fftsize/2+7]=1
  elif spotpsfdim==4:
    spotpsf=Numeric.zeros((nsubx,nsubx,fftsize,fftsize),"f")
    spotpsf[:,:,fftsize/2-5:fftsize/2+5,fftsize/2-5:fftsize/2+5]=1
    spotpsf[:,:,fftsize/2-7:fftsize/2-3,fftsize/2-7:fftsize/2-3]=1
    spotpsf[:,:,fftsize/2+3:fftsize/2+7,fftsize/2+3:fftsize/2+7]=1
    spotpsf[2,2,fftsize/2-4:fftsize/2+4,fftsize/2-4:fftsize/2+4]=0
    #spotpsf[2,2]=0.
  else:
    spotpsf=None
  if type(sig)==Numeric.ArrayType:
    sig=sig.flat
  bimg=Numeric.zeros((nsubx,nsubx,nimg,nimg),Numeric.Float32)
  print "Initialising"
  cc=util.centcmod.centcmod(nthreads,nsubx,ncen,fftsize,nimg,phasesize,readnoise,readbg,addPoisson,noiseFloor,sig,skybrightness,calsource,pxlPower,nintegrations,seed,reorderedPhs,pup,spotpsf,cents,bimg)
  print "Running"
  t=cc.run(calsource)
  print "Time taken",t

  def makeImage(bimg):
    img=Numeric.zeros((nsubx*nimg,nsubx*nimg),Numeric.Float32)
    for i in range(nsubx):
      for j in range(nsubx):
        img[i*nimg:(i+1)*nimg,j*nimg:(j+1)*nimg]=bimg[i,j]
    return img

  gist.fma();gist.pli(makeImage(bimg))
