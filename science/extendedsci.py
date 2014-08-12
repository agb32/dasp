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


class extsci(base.aobase.aobase):
    def __init__(self,parentDict,config,args={},forGUISetup=0,debug=None,idstr=None):
        base.aobase.aobase.__init__(self,parentDict,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.parentDict=parentDict
        self.gridx=self.config.getVal("shsGridX")
        self.gridy=self.config.getVal("shsGridY")
        if len(self.parentDict.keys())!=self.gridx*self.gridy:
            raise Exception("Wrong number of parents...")
        self.refImage=self.config.getVal("shsRefImage")#The reference image to be used when convolving the PSFs together.
        self.nimg=self.config.getVal("scinimg")
        self.nsection=self.config.getVal("n_section")
        self.insimg=numpy.empty([self.nimg,self.nimg])


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
                self.processData()
        else:
            self.dataValid=0

                
    def processData(self):
        """The point SHS images are now stored in self.inputData[(y,x)].
        So, combine these together in some way, to give the extended SHS image.
        Probably a shift and sum (scaled by the reference image), or maybe just a sum (scaled by the reference image).
        """
        #print "todo"

        sci_nimg=self.config.getVal("sci_nimg")
        finalimg=self.config.getVal("finalimg")
        sn_section=self.config.getVal("sn_section")
        ss_ovlp=self.config.getVal("ss_ovlp")
        tmpres=numpy.empty([finalimg,finalimg])
        tmpres[:]=0
        tmppsf=numpy.empty([sci_nimg+sn_section+ss_ovlp*2-1,sci_nimg+sn_section+ss_ovlp*2-1])
        tmppsf[:]=0
        ssfc=self.config.getVal("ssfcforsci")
        tover=self.config.getVal("scit_over")
        n_gap=(sci_nimg)/2
        ctv=tmpres*0
        
        for y in range(self.gridy):
            for x in range(self.gridx):
                fftss=ssfc[x+self.gridy*y].copy()
                tmpdata=self.parentDict["%d"%(x+y*self.gridx)].insimg
                gap=(tmpdata.shape[0]-sci_nimg-1)/2
               # print "gap:%d"%gap
                tmppsf[0:sci_nimg,0:sci_nimg]=tmpdata[gap:gap+sci_nimg,gap:gap+sci_nimg]
                fftc=FFT.fft2(tmppsf)*fftss
                conv=FFT.ifft2(fftc).real
                
                boxy0=max(y*sn_section-tover,0)
                boxy1=min((y+1)*sn_section+tover,finalimg)
                boxx0=max(x*sn_section-tover,0)
                boxx1=min((x+1)*sn_section+tover,finalimg)
                #print "%d,%d,%d,%d"%(boxy0,boxy1,boxx0,boxx1)
                
                imboxy0=boxy0-(y*sn_section-tover)+n_gap+ss_ovlp-tover
                imboxy1=imboxy0+(boxy1-boxy0)
                imboxx0=boxx0-(x*sn_section-tover)+n_gap+ss_ovlp-tover
                imboxx1=imboxx0+(boxx1-boxx0)
                tmpres[boxy0:boxy1,boxx0:boxx1]+=conv[imboxy0:imboxy1,imboxx0:imboxx1]
                ctv[boxy0:boxy1,boxx0:boxx1]+=1
        self.insimg=tmpres/ctv
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
        txt="""<plot title="extendedshs Output data%s" cmd="data=%s.insimg" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        return txt
	
