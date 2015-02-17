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
