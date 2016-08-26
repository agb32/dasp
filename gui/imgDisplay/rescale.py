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
import types,numpy
class rescale:
    """A class to rescale, resize, reposition an image, ready for plotting"""
    def __init__(self,size=None,zoom=[1,1]):
        """size is a 2 tuple, the x and y dimensions"""
        
        if size!=None:
            self.output=numpy.zeros(size,numpy.uint8)
        else:
            self.output=None


    def rescale(self,img,size=None,pos=None,zoom=None,minval=0,maxval=255,auto=None,divFactor=1.):
        """size is a 2-tuple
        pos is a 2-tuple (x and y position in img at which to start
        zoom is a 2-tuple (x and y zoom factors)
        if auto: scale image using full range of 8 bits.
        Otherwise use minval,maxval after dividing the image by divfactor
        This should probably be speeded up in c eventually.
        """
        if size==None:
            if type(self.output)!=types.NoneType:
                size=self.output.shape
            else:
                size=img.shape
        if type(self.output)!=types.NoneType:
            if self.output.shape!=size:
                self.output=numpy.zeros(size,numpy.uint8)
            else:
                self.output[:]=0
        else:
            self.output=numpy.zeros(size,numpy.uint8)
        if type(zoom)==types.NoneType:
            zoom=self.zoom
        if type(pos)==types.NoneType:
            pos=self.pos
        input=numpy.zeros(((size[0]+zoom[0]-1)/zoom[0],(size[1]+zoom[1]-1)/zoom[1]),img.dtype)
        tmp=img[pos[0]:pos[0]+size[0]/zoom[0],pos[1]:pos[1]+size[1]/zoom[1]]
        input[0:tmp.shape[0],0:tmp.shape[1]]=tmp[:,:]
        inputmax=max(input.flat)
        inputmin=min(input.flat)
        inputdiff=float(inputmax-inputmin)
        if inputdiff==0.0:
            inputdiff=1.0
        if auto:
            #print size,self.output.shape,input.shape
            input-=inputmin
            inputdiff/=255.
            input/=inputdiff
            for i in range(size[0]):
                for j in range(size[1]):
                    #print i,j,i/zoom[0],j/zoom[1]
                    self.output[i,j]=input[i/zoom[0],j/zoom[0]]
        else:
            divFactor/=255.
            input/=divFactor
            input=numpy.where(input>maxval,maxval,input)
            input=numpy.where(input<minval,minval,input)
            for i in range(size[0]):
                for j in range(size[1]):
                    self.output[i,j]=input[i/zoom[0],j/zoom[0]]
        return self.output


    def minmaxNum(self,typecode):
        """minimum and maximum possible values for integer datatypes"""
        rtval=None
        if typecode=='b':
            rtval=(0,255)
        elif typecode=='1':
            rtval=(-128,127)
        elif typecode=="s":
            rtval=(-32768,32767)
        elif typecode=="w":
            rtval=(0,65536)
        elif typecode=="u":
            rtval=(0,4294967295)
        elif typecode=="i":
            rtval=(-2147483648,2147483647)
        return rtval
