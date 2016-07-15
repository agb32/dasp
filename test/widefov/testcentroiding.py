import numpy
import util.centroid
import util.tel
nsubx=10
npup=80
nimg=8
image=numpy.zeros((80,80),"f")
#image[4::8,4::8]=1

fftsize=16
phasesize=8
ncen=8
pup=util.tel.Pupil(npup,npup/2,0,nsubx)
c=util.centroid.centroid(nsubx,pup,fftsize=fftsize,binfactor=None,phasesize=phasesize,nimg=nimg,ncen=ncen,inputImage=image)
c.easy()
c.runSlopeCalc({"cal_source":0})

#To change the ref slopes to a given reference image, set cameraInput equal to the image (if 0, will be a 0 image).
c.centcmod.update(util.centcmod.REFCENTS,None)
c.takeReference({"cal_source":1},cameraInput=0)
