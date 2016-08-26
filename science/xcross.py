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
"""A modular to calculate the corss-correlation for extended WFS
Base on Lofdahl A&A 2010 and the code of Dr. ANC Kellerer

Author: G. Zhao, NIAOT

Version 1.0, 24 Jan 2014
 Just a normal one...
 things to do: C code?
 
"""


import util.FITS
import matplotlib.pyplot as plt
import numpy.fft as fft
import numpy
import string
import base.aobase
import scipy.ndimage as ndimg
import math

class xcross(base.aobase.aobase):
    def __init__(self,parentDict,config,args={},forGUISetup=0,debug=None,idstr=None):
        base.aobase.aobase.__init__(self,parentDict,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.parentDict=parentDict
        nimg=self.config.getVal("n_wfs_subapt")
        self.nimg=nimg
        self.nsubx=self.config.getVal("wfs_nsubx")
        self.nsuby=self.nsubx
        self.errorall=0
        self.shiftall=0
        
        # read the correlation algorithms mode
        # CA_mode need to be in the format of "SDF0"
        # the first three letters are the name of CA, which can be "SDF, CFI, CFF, ADF, ADS, NCF"
        # see Table 1 of Lofdahl 2010 for the detail
        # NCF computes the normalized cross-correlation, same as the Matlab function normxcorr2()
        # the last number of CA_mode is subpixel interpolation methos
        # 0,1,2,3 means 2QI, 2LS, 1QI, 1LS
        
        mode0=self.config.getVal("CA_mode", default="""NCF1""")
        mode=mode0[0:3]
        
        if mode == "SDF":
            self.init_mode = 0
            self.cor_mode= 0
        elif mode == "CFI":
            self.init_mode = 1
            self.cor_mode= 2
        elif mode == "CFF":
            self.init_mode = 1
            self.cor_mode= 1
        elif mode =="ADF":
            self.init_mode = 0
            self.cor_mode= 3
        elif mode == "ADS":
            self.init_mode = 0
            self.cor_mode= 4
        elif mode == "NCF":
            self.init_mode = 0
            self.cor_mode= 5
        else:
            # or U can use something like "M002"
            self.init_mode = ord(mode0[1])-ord('0')
            self.cor_mode = ord(mode0[2])-ord('0')

        self.intp_mode = ord(mode0[3])-ord('0')
        
        #A number of iteration, in some case, the reference image need to be update every several iteration

        self.iter=0
        
        ny=nimg
        nx=nimg
        self.sz_img=[ny,nx]
        # the reference image, use self.load_ref to reload one
        
        #some globle constants used in every iteration 
        self.linex=numpy.linspace(0,nx-1,nx)
        self.liney=numpy.linspace(0,ny-1,ny)
        w1x=numpy.zeros([1,nx])
        w1y=numpy.zeros([ny,1])
        a=0.53836
        w1x[0,:]=a+(a-1)*numpy.cos(2.*numpy.pi*self.linex/(nx-1))
        w1y[:,0]=a+(a-1)*numpy.cos(2.*numpy.pi*self.liney/(ny-1))
        self.wdw=numpy.dot(w1y,w1x)
        w1x[0,:]=self.linex
        w1y[:,0]=self.liney
        self.xgrid=numpy.dot((w1y*0+1),w1x)
        self.ygrid=numpy.dot(w1y,(w1x*0+1))
        
    def generateNext(self):
        """Compression main loop - compress along a line of sight to get WFS images values."""
        if self.generate==1:
            if self.newDataWaiting:
                nin=0
                if self.parentDict.dataValid==1:
                    nin+=1
                else:
                    print "Xcross: Waiting for data from wfs, but not valid"
                if nin>0:
                    self.dataValid=1
                else:
                    self.dataValid=0
            if self.dataValid:
                self.processData()
        else:
            self.dataValid=0
            
    def processData(self):
        inputimg=self.parentDict.shimg
        nimg=self.nimg
        
        if self.iter == 0:
            # first time run the modular
            # which one to do the modular

            self.ifimage=self.config.getVal("imagebod",default=1)
            self.full=[]
            maxsub=[]
            for pupx in range(self.nsubx):
                for pupy in range(self.nsuby):
                    tmpdata=inputimg[pupy*nimg:(pupy+1)*nimg,pupx*nimg:(pupx+1)*nimg]

                    sumtmp=tmpdata.sum()
                    if sumtmp > self.ifimage:
                        self.full.append([pupy,pupx])
                        maxsub.append(sumtmp)
            ref_pos=self.config.getVal("ref_pos",default=[-1,-1])
            if ref_pos[0] < 0:
                imax=maxsub.index(max(maxsub))
                ref_pos=self.full[imax]
                
            ref=inputimg[ref_pos[0]*nimg:(ref_pos[0]+1)*nimg,ref_pos[1]*nimg:(ref_pos[1]+1)*nimg]
            self.ref=self.load_ref(ref)
            self.ref_shift=[0,0]
            self.slope_shift=self.parentDict.outputData[ref_pos[0],ref_pos[1],:]
        self.iter +=1
        
        allslope=numpy.zeros([self.nsuby,self.nsubx,2])
        for i in range(len(self.full)):
            pos=self.full[i]
            img=inputimg[pos[0]*nimg:(pos[0]+1)*nimg,pos[1]*nimg:(pos[1]+1)*nimg]
            s1=self.f_subpixel(self.initial_img(img))

            allslope[pos[0],pos[1],:]=s1+self.slope_shift

        #print allslope[:,:,0]
        #print self.parentDict.outputData[:,:,0]
        shwow=numpy.zeros([self.nsuby,self.nsubx*2])
        shwow[:,0:self.nsubx]=allslope[:,:,0]
        shwow[:,self.nsubx:self.nsubx*2]=self.parentDict.outputData[:,:,0]
        plt.imshow(shwow)
        plt.show()
        error=((allslope[:,:,:]-self.parentDict.outputData[:,:,:])**2).sum()/(37*2.)
        shift=(self.parentDict.outputData[:,:,:]**2).sum()/(37*2.)
        self.errorall=self.errorall+error
        self.shiftall+=shift
        if self.iter == 100:
            print "error=%f, avg_shift=%f"%(math.sqrt(self.errorall/100), math.sqrt(self.shiftall/100))
       # raw_input()
        
    def avg(self,array):
        return array.sum()/array.size
        
    def f_subpixel(self,img):

        nx=self.sz_img[1]
        ny=self.sz_img[0]

        amp=int(min(nx/4.,ny/4.,10.))
        if self.cor_mode != 1:
            bxsz=int(min(nx/2,ny/2))
            i0=int((nx-bxsz)/2.);
            i1=i0+bxsz-1;
            j0=int((ny-bxsz)/2.);
            j1=j0+bxsz-1;
            sm_ref=self.ref[i0:i1,j0:j1]
            
            cor=numpy.zeros([2*amp+1,2*amp+1])
            
            for di in range(-amp,amp+1):
                for dj in range(-amp,amp+1):
                    im=img[i0+di:i1+di,j0+dj:j1+dj]
                    if self.cor_mode == 0:
                        cor[di+amp,dj+amp]=((im-sm_ref)**2).sum()
                        #if di==0 and dj ==0:
                            #print "%f"%((im-sm_ref)**2).sum()
                            
                    elif self.cor_mode == 2:
                        cor[di+amp,dj+amp]=-(im*sm_ref).sum()
                    elif self.cor_mode == 3 or self.cor_mode == 4:
                        cor[di+amp,dj+amp]=(abs(im-sm_ref)).sum()
                    elif self.cor_mode == 5:
                        down=((im-self.avg(im))**2).sum()*((sm_ref-self.avg(sm_ref))**2).sum()
                        cor[di+amp,dj+amp]=-((im-self.avg(im))*(sm_ref-self.avg(sm_ref))).sum()/math.sqrt(down)

            s=cor.argmin()
            y0=int(s/cor.shape[1])
            x0=s%cor.shape[1]

            y1=y0-amp
            x1=x0-amp

        elif self.cor_mode == 1:

            cor=fft.ifft2(fft.fft2(img*self.wdw) * self.ref)
            #cor=fft.ifft2(numpy.conjugate(self.ref) * self.ref)

            cor=-abs(cor)**2

            s=cor.argmin()
            y0=int(s/cor.shape[1])
            x0=s%cor.shape[1]

            if x0<nx/2:
                x1=x0
            else:
                x1=-nx+x0

            if y0<ny/2:
                y1=y0
            else:
                y1=-ny+y0

        cc=numpy.zeros([3,3])
        cx=cor.shape[1]
        cy=cor.shape[0]

        for i in range(3):
            for j in range(3):
                cc[j,i]=cor[(y0-1+j)%cy,(x0-1+i)%cx]
        
        if self.cor_mode == 4:
            cc=cc**2
        
        if self.intp_mode == 0 or self.intp_mode ==2:

            a2=(sum(cc[2,:])-sum(cc[0,:]))/6.
            a3=(sum(cc[2,:])-2*sum(cc[1,:])+sum(cc[0,:]))/6.
            a4=(sum(cc[:,2])-sum(cc[:,0]))/6.
            a5=(sum(cc[:,2])-2*sum(cc[:,1])+sum(cc[:,0]))/6.
            a6=(cc[2,2]-cc[0,2]-cc[2,0]+cc[0,0])/4.

        elif self.intp_mode == 1 or self.intp_mode ==3:

            a2=(cc[2,1]-cc[0,1])/2.
            a3=(cc[2,1]-2*cc[1,1]+cc[0,1])/2.
            a4=(cc[1,2]-cc[1,0])/2.
            a5=(cc[1,2]-2*cc[1,1]+cc[1,0])/2.
            a6=(cc[2,2]-cc[0,2]-cc[2,0]+cc[0,0])/4.
            
        if self.intp_mode == 0 or self.intp_mode == 1: 
            a6635=(a6*a6-4.*a3*a5)
            addy=(2.*a2*a5-a4*a6)/a6635
            addx=(2.*a3*a4-a2*a6)/a6635
            
            if abs(addy) <= 1:
                y2=y1+addy
            else:
                y2=y1-a2/a3/2.
            
            if abs(addx) <= 1:
                x2=x1+addx
            else:
                x2=x1-a4/a5/2.
        else:
            y2=y1-a2/a3/2.
            x2=x1-a4/a5/2.
        
        return [x2,y2]
        
    def initial_img(self,ref):
        ref1=ref/self.avg(ref)
        if self.init_mode == 0:
            ref1-=self.avg(ref1)
        elif self.init_mode == 1:
            bx=numpy.polyfit(self.linex,ref1.sum(axis=0)/self.sz_img[0],1)
            ref1 -= bx[1]+bx[0]*self.xgrid
            by=numpy.polyfit(self.liney,ref1.sum(axis=1)/self.sz_img[1],1)
            ref1 -= by[1]+by[0]*self.ygrid
        return ref1
        
    def load_ref(self,ref):
        ref1=self.initial_img(ref)
        if self.cor_mode == 1:
            ref1=numpy.conjugate(fft.fft2(ref1*self.wdw))
        return ref1
        
    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr[0]==None or self.idstr[0]=="":
            id=""
        #else:
        #    id=" (%s)"%self.idstr[0]
        #txt="""<plot title="extendedshs Output data%s" cmd="data=%s.shimg" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        return ""
