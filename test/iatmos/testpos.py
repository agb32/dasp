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
#Tries to work out rotations etc so that guide stars maintain the correct positional locations with each other.
import numpy
import cmod.iscrn
import pylab
fsize=20
rdata=numpy.random.random((fsize,fsize))
from scipy.interpolate import interp2d
interp=interp2d(numpy.arange(fsize),numpy.arange(fsize),rdata,kind="cubic")
scrn=interp(numpy.arange(1000)/999.*(rdata.shape[0]-1),numpy.arange(1000)/999.*(rdata.shape[1]-1)).copy()
out=numpy.zeros((500,500),"f")
dirn=45.
istr=cmod.iscrn.initialiseInterp(scrn,None,-dirn,out,1.,0,1,1)

pltpos=[8,3,2,5,10,11]
outlist=[]
for i in range(6):
    rotangle=dirn+180.
    xpos=100*numpy.cos((i*60+rotangle)*numpy.pi/180.)
    ypos=100*numpy.sin((i*60+rotangle)*numpy.pi/180.)
    out[:]=0
    cmod.iscrn.rotShiftWrapSplineImageThreaded(istr,xpos,ypos,0)
    outlist.append(out.copy())
    pylab.subplot(3,4,pltpos[i])
    pylab.imshow(outlist[-1])
pylab.show()
