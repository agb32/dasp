"""A module to created extended object WFS images


"""

#from science.xinterp_dm import interpolateSpline,interpolateBicubic
import numpy
import time
import base.aobase
import util.FITS
import scipy
import scipy.fftpack as FFT
import matplotlib.pyplot as plt
from scipy.ndimage import interpolation


class SHS(base.aobase.aobase):
    def __init__(self,parentDict,config,args={},forGUISetup=0,debug=None,idstr=None):
        base.aobase.aobase.__init__(self,parentDict,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.parentDict=parentDict
        self.dmObj=self.config.getVal("dmObj")
        self.gridx=self.config.getVal("shsGridX")
        self.gridy=self.config.getVal("shsGridY")
        if len(self.parentDict.keys())!=self.gridx*self.gridy:
            raise Exception("Wrong number of parents...")
        self.refImage=self.config.getVal("shsRefImage")#The reference image to be used when convolving the PSFs together.
        self.nimg=self.config.getVal("wfs_nimg")
        self.nsubx=self.config.getVal("wfs_nsubx")
        self.nsection=self.config.getVal("n_section")
        self.shimg=numpy.zeros((self.nsection*self.nsubx*self.gridx,self.nsection*self.nsubx*self.gridy),numpy.float32)
        self.ovlp=self.config.getVal("s_ovlp")
        self.outputData=numpy.empty([self.nsubx,self.nsubx,2])
        self.tovlp=self.config.getVal("t_ovlp")
        self.loop=0
        self.outputkey=0


    def generateNext(self):
        """Compression main loop - compress along a line of sight to get DM actuator values."""
        if self.generate==1:
            if self.newDataWaiting:
                nin=0
                for key in self.parentDict.keys():
                    if self.parentDict[key].dataValid==1:
                        nin+=1
                    else:
                        print "tomoRecon: Waiting for data from wfs, but not valid"
                if nin>0:
                    if nin==len(self.parentDict.keys()):
                        self.dataValid=1
                    else:
                        print "tomoRecon: Warning - got some data but not all, setting dataValid=0"
                        self.dataValid=0
                else:
                    self.dataValid=0
            if self.dataValid:
                self.getInput()
                self.processData()
        else:
            self.dataValid=0


    def getInput(self):
        self.inputData={}
        for y in range(self.gridy):
            for x in range(self.gridx):
               # print type(self.parentDict["%d"%(x+y*self.gridx)].subimg)
                #print self.parentDict["%d"%(x+y*self.gridx)].subimg.shape
    #            print self.parentDict["%d"%(x+y*self.gridx)].shimg2.shape
    #            print self.parentDict["%d"%(x+y*self.gridx)].outputData.shape
                self.inputData[(y,x)]=self.parentDict["%d"%(x+y*self.gridx)].outputData
    #            plt.imshow(self.parentDict["%d"%(x+y*self.gridx)].shimg2)
    #            plt.show()
    #            print self.parentDict["%d"%(x+y*self.gridx)].outputData.shape
                
    def processData(self):
        """The point SHS images are now stored in self.inputData[(y,x)].
        So, combine these together in some way, to give the extended SHS image.
        Probably a shift and sum (scaled by the reference image), or maybe just a sum (scaled by the reference image).
        """
        #print "todo"
        self.shimg[:]=0
        nsection=self.nsection
        nimg=self.nimg
        nbigimg=self.nsection*self.gridx
        nallimg=nbigimg*self.nsubx
        tmpout=numpy.empty([nbigimg*self.nsubx,nbigimg*self.nsubx])
        tmpout[:]=0
        ctv=tmpout*0
        tover=self.tovlp
            
        ssfc=self.config.getVal("ssfc")
        #tmpss=numpy.empty([nimg+nsection-1,nimg+nsection-1])
        #tmppsf=tmpss.copy()
        tmppsf=numpy.empty([nimg+nsection+self.ovlp*2-1,nimg+nsection+self.ovlp*2-1])
        shiftpsf=tmppsf.copy()
        neach=nimg+nsection+self.ovlp*2-1
        n_gap=(nimg-1)/2
        self.outputData=0
        for y in range(self.gridy):
            for x in range(self.gridx):
                
                fftss=ssfc[x+self.gridy*y].copy()
                tmpdata=self.parentDict["%d"%(x+y*self.gridx)].drawCents(0)#shimg2.copy()
                #shiftvalue=self.parentDict["%d"%(x+y*self.gridx)].outputData.copy()
                self.outputData+=self.parentDict["%d"%(x+y*self.gridx)].outputData
                
                    
                #print shiftvalue
                for pupx in range(self.nsubx):
                    for pupy in range(self.nsubx):
                        tmppsf[:]=0
                        tmppsf[0:nimg,0:nimg]=tmpdata[pupy*nimg:(pupy+1)*nimg,pupx*nimg:(pupx+1)*nimg]
                        fftc=FFT.fft2(tmppsf)*fftss
                        conv=FFT.ifft2(fftc).real
                        boxy0=max(y*nsection+pupy*nbigimg-tover,nbigimg*pupy)
                        boxy1=min((y+1)*nsection+pupy*nbigimg+tover,nbigimg*(pupy+1))
                        boxx0=max(x*nsection+pupx*nbigimg-tover,nbigimg*pupx)
                        boxx1=min((x+1)*nsection+pupx*nbigimg+tover,nbigimg*(pupx+1))
                        #print "%d,%d,%d,%d"%(boxy0,boxy1,boxx0,boxx1)
                        
                        imbox=n_gap+self.ovlp
                        imboxy0=boxy0-(y*nsection+pupy*nbigimg-tover)+n_gap+self.ovlp-tover
                        imboxy1=imboxy0+(boxy1-boxy0)
                        imboxx0=boxx0-(x*nsection+pupx*nbigimg-tover)+n_gap+self.ovlp-tover
                        imboxx1=imboxx0+(boxx1-boxx0)
                        #print "%d,%d,%d,%d"%(imboxy0,imboxy1,imboxx0,imboxx1)
                        #raw_input()
                        tmpout[boxy0:boxy1,boxx0:boxx1]+=conv[imboxy0:imboxy1,imboxx0:imboxx1]
                        ctv[boxy0:boxy1,boxx0:boxx1]+=1
                        #for ii in range(nsection):#+self.ovlp*2):
                            #yinbig=y*nsection+pupy*nbigimg+ii#-self.ovlp
                            #for jj in range(nsection):#+self.ovlp*2):
                                #xinbig=x*nsection+pupx*nbigimg+jj#-self.ovlp
                                #if yinbig >= pupy*nbigimg and yinbig < (pupy+1)*nbigimg and xinbig >= pupx*nbigimg and xinbig < (pupx+1)*nbigimg:
                                #tmpout[yinbig,xinbig]+=conv[ii+n_gap+self.ovlp,jj+n_gap+self.ovlp]
                #tmpout1=FFT.fftshift(FFT.fft2(tmpout))
        tmpout1=tmpout/ctv
        self.loop=self.loop+1
        Nframe=self.config.getVal("Nframe_of_outputfile",default=-1)
        if self.loop <= Nframe and Nframe > 0:
            file_name_of_fits=self.config.getVal("Name_of_output_fits_file")
            util.FITS.WriteCube(tmpout1,file_name_of_fits,Nframe)
            util.FITS.WriteCube(tmpdata,"%spsf"%file_name_of_fits,Nframe)
            print self.loop
            if self.loop == Nframe:
                print "Finish output extendedSHS output fits file for %d frames"%Nframe
        self.outputData=self.outputData/(self.gridx*self.gridy)
        #print self.outputData.shape
        """
        self.outputkey=self.outputkey+1
        print "outputkey=%d"%self.outputkey
        if self.outputkey == 1:
            self.filek=open('/Users/gangzhao/ProgramFiles/DASP/solar_aglae/slope.txt','w')
        if self.outputkey < 5001:
            tmp=self.outputData[0]
            print type(tmp)
            print tmp.shape
            self.filek.write("%f,%f\n"%(tmp[0,0],tmp[0,1]))
            print "finish"
        if self.outputkey == 5001:
            self.filek.close()
       # for ii in range(self.nsubx):
        #    for x in range(self.nsubx):
       """  
                
                         
        self.shimg[:]=tmpout1/1e7
        """or key in self.parentDict.keys():
        self.outputData[:]+=self.parentDict[key].outputData
        """
        
        
    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr[0]==None or self.idstr[0]=="":
            id=""
        else:
            id=" (%s)"%self.idstr[0]
        txt="""<plot title="extendedshs Output data%s" cmd="data=%s.shimg" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        return txt
	
