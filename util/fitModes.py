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
"""The methods implemented here can be used to fit mirror shapes to the mirror,
ie compute the actuator values needed to best fit."""
import numpy as na
from util.dm import MirrorSurface#interpolateBicubic,interpolateSpline

class fitModes:
    def __init__(self,npup,nact,nmodes,modes,interpType,actCoupling,actFlattening):
        #self.config=config
        self.interpType=interpType
        self.actCoupling=actCoupling
        self.actFlattening=actFlattening
        self.modes=na.array(modes).astype(na.float32)
        self.nmodes=nmodes
        self.npup=npup
        #self.dm=science.xinterp_dm.dm(None,config)
        self.nact=nact

        self.actarr=na.zeros((self.nmodes,self.nact,self.nact),na.float32)
        self.interpolatedAct=na.zeros((self.nmodes,self.npup,self.npup),na.float32)
        self.mirrorSurface=MirrorSurface(self.interpType,self.npup,self.nact,actCoupling=actCoupling,actFlattening=actFlattening)
    def interpolate(self,actmap):
        dmphs=self.mirrorSurface.fit(actmap)#interp(actmap,self.interpType,self.npup,self.actCoupling,self.actFlattening,self.nact)
##         if self.interpType=="bicubic":
##             dmphs=interpolateBicubic(actmap,self.npup,self.actCoupling,self.actFlattening)
##         else:#bicubic spline
##             dmphs=interpolateSpline(actmap,self.npup,self.nact)
        return dmphs

    def fit(self,actarr=None):
        """The modes is a 3d array, (nmodes, npup, npup), the modes what are to
        be fitted onto the mirror.
        Note, the modes should not be constrained by a pupil as this will allow
        better fitting (I think).
        At the moment, doesn't do any simplexing or anything like that...
        """
        if type(actarr)==type(None):
            step=float(self.npup-1)/float(self.nact-1)
            for i in range(self.nmodes):
                for y in range(self.nact):
                    for x in range(self.nact):
                        xx=int(x*step+0.5)
                        yy=int(y*step+0.5)
                        self.actarr[i,y,x]=self.modes[i,yy,xx]
            actarr=self.actarr
        for i in range(self.nmodes):
            self.interpolatedAct[i]=self.interpolate(actarr[i])
def interp(actmap,interpType,npup,actCoupling,actFlattening,nact):
    m=MirrorSurface(interpType,npup,actmap.shape[0],actCoupling=actCoupling,actFlattening=actFlattening)
    dmphs=m.fit(actmap)
    #if interpType=="bicubic":
    #    dmphs=interpolateBicubic(actmap,npup,actCoupling,actFlattening)
    #else:#bicubic spline
    #    dmphs=interpolateSpline(actmap,npup)
    return dmphs
    
def testLinearity(actarr1,actarr2,interpType="bicubic",npup=64,nact=9,actCoupling=0.1,actFlattening=0.25):
    """tests the linearity between 2 modes"""
    if interpType=="bicubic":
        phs1=interpolateBicubic(actarr1,npup,actCoupling,actFlattening)
        phs2=interpolateBicubic(actarr2,npup,actCoupling,actFlattening)
        phs3=interpolateBicubic(actarr1+actarr2,npup,actCoupling,actFlattening)
    else:
        phs1=interpolateSpline(actarr1,npup,nact)
        phs2=interpolateSpline(actarr2,npup,nact)
        phs3=interpolateSpline(actarr1+actarr2,npup,nact)
    phs4=phs1+phs2
    phs5=phs3-phs4
    return phs3,phs4,phs5
