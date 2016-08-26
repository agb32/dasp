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
