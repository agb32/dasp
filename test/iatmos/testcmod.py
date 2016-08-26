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
import cmod.iscrn
from scipy.interpolate import interp2d

import numpy

import time

import pylab

def rotateShiftWrapSplineImage(img,dimg,deg,shiftx,shifty,wrappoint,out,r=1.):
    """Rotates an image by deg degrees, shifts centre by x,y and puts result into out, which can be different size from img.
    wrappoint is the point at which img is wrapped, being the oldest phase.
    dimg is the gradient of img in y direction, ie approx (img[2:]-img[:-2])/2.   (with the edge cases included too).
    r is a zoom factor - e.g 0.5 gives zoom of 2.
"""
    ang=deg*numpy.pi/180.
    s=r*numpy.sin(ang)
    c=r*numpy.cos(ang)
    dim=out.shape
    imgdim=img.shape
    sx=shiftx
    sy=shifty
    maxx=0
    maxy=0
    minx=None
    miny=None
    points=numpy.zeros((4,),"f")
    outofrange=numpy.zeros((4,),"i")
    for yy in range(dim[0]):
        y=yy-dim[0]/2.+.5
        for xx in range(dim[1]):
            x=xx-dim[1]/2.+.5
            yold=-s*x+c*y+imgdim[0]/2.-.5-sy+wrappoint
            xold=c*x+s*y+imgdim[1]/2.-.5-sx
            x1=int(numpy.floor(xold))
            x2=x1+1
            x0=x1-1
            x3=x1+2
            #First, we need to compute 4 splines in the y direction.  These are then used to compute in the x direction, giving the final value.
            y1=int(numpy.floor(yold))
            if y1>=imgdim[0]:#wrap it
                y1-=imgdim[0]
                yold-=imgdim[0]
            y2=y1+1
            xm=xold-x1
            ym=yold-y1
            if y2==imgdim[0]:
                y2=0

            #X1=0, X2=1.  
            #t=ym
            for i in range(4):#at x1-1, x1, x2, x2+1.
                if x1+i-1>=0 and x1+i-1<img.shape[1]:
                    k1=dimg[y1,x1+i-1]
                    k2=dimg[y2,x1+i-1]
                    Y1=img[y1,x1+i-1]
                    Y2=img[y2,x1+i-1]
                    a=k1-(Y2-Y1)#k1*(X2-X1)-(Y2-Y1)
                    b=-k2+(Y2-Y1)#-k2*(X2-X1)+(Y2-Y1)
                    points[i]=(1-ym)*Y1+ym*Y2+ym*(1-ym)*(a*(1-ym)+b*ym)
                    outofrange[i]=0
                else:
                    outofrange[i]=1
            #and now interpolate in X direction (using points).
            if outofrange[0]:
                k1=points[2]-points[1]
            else:
                k1=(points[2]-points[0])*.5
            if outofrange[3]:
                k2=points[2]-points[1]
            else:
                k2=(points[3]-points[1])*.5
            if outofrange[1] or outofrange[2]:
                raise Exception("Out of range")
            #t=xm.
            #X1=0, X2=1
            Y1=points[1]
            Y2=points[2]
            a=k1-(Y2-Y1)#k1*(X2-X1)-(Y2-Y1)
            b=-k2+(Y2-Y1)#-k2*(X2-X1)+(Y2-Y1)
            val=(1-xm)*Y1+xm*Y2+xm*(1-xm)*(a*(1-xm)+b*xm)
            
            out[yy,xx]=val
    return out




size=1280
ssize=int(numpy.ceil(1280*numpy.sqrt(2)+1))

rnd=numpy.random.random((ssize//100,ssize//100))
interp=interp2d(numpy.arange(ssize//100),numpy.arange(ssize//100),rnd,kind="cubic")
rndinterp=interp(numpy.arange(ssize)/(ssize-1.)*(rnd.shape[1]-1.),numpy.arange(ssize)/(ssize-1.)*(rnd.shape[0]-1.))
print rndinterp.shape,ssize


screen=rndinterp.copy()#numpy.random.random((ssize,ssize)).astype("d")
dscreen=numpy.empty((ssize,ssize),"d")
dscreen[1:-1]=(screen[2:]-screen[:-2])*0.5
dscreen[0]=screen[1]-screen[0]
dscreen[-1]=screen[-1]-screen[-2]

out=numpy.zeros((size,size),"f")
out2=numpy.zeros((size,size),"f")
m1=m2=m3=m4=(10,0,0,0)
for nthr in range(4,5):
    for nx in range(1,9):
        for ny in range(1,9):
            interpStruct=cmod.iscrn.initialiseInterp(screen,dscreen,80.,out,1.,nthr,nx,ny)
            interpStructNoGrad=cmod.iscrn.initialiseInterp(screen,None,80.,out2,1.,nthr,nx,ny)


            t1=time.time()
            cmod.iscrn.rotShiftWrapSplineImageThreaded(interpStruct,0.,0.,0)
            #rotateShiftWrapSplineImage(screen,dscreen,45.,0.,0.,0,out,1.)
            t2=time.time()
            cmod.iscrn.rotShiftWrapSplineImageThreaded(interpStructNoGrad,0.,0.,0)
            t3=time.time()
            cmod.iscrn.rotShiftWrapSplineImageThreaded(interpStruct,0.,0.,0)
            t4=time.time()
            cmod.iscrn.rotShiftWrapSplineImageThreaded(interpStructNoGrad,0.,0.,0)
            t5=time.time()
            if t2-t1<m1[0]:
                m1=(t2-t1,nthr,nx,ny)
            if t3-t2<m2[0]:
                m2=(t3-t2,nthr,nx,ny)
            if t4-t3<m3[0]:
                m3=(t4-t3,nthr,nx,ny)
            if t5-t4<m4[0]:
                m4=(t5-t4,nthr,nx,ny)

            print "%d %d %d Took %g %g %g %g"%(nthr,nx,ny,t2-t1,t3-t2,t4-t3,t5-t4)
print numpy.min(screen),numpy.max(screen),screen.mean(),numpy.min(out),numpy.max(out),out.mean()
print m1,m2,m3,m4
pylab.subplot(2,2,1)#figure(1)
pylab.imshow(screen)
pylab.gca().invert_yaxis()
pylab.subplot(2,2,2)
#pylab.figure(2)
pylab.imshow(out)
pylab.gca().invert_yaxis()
pylab.subplot(2,2,3)
pylab.imshow(out2)
pylab.gca().invert_yaxis()
pylab.subplot(2,2,4)
pylab.imshow(out2-out)
pylab.gca().invert_yaxis()
pylab.show()




(0.027168989181518555, 4, 7, 4)
(0.03361701965332031, 4, 6, 2)
(0.02780890464782715, 4, 8, 1)
(0.03362703323364258, 4, 8, 3)
