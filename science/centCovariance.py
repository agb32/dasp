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
import util.FITS
import base.aobase
class cov(base.aobase.aobase):
    """
    Module to compute centroid covariance.  Needs 2 parents, both of which return centroids, one with no shot or read noise, and the other with noise.  When computing covariances (actually, variances), the DM should be flat, but the atmosphere should be running.  The variances of each centroider are computed (over lots of iterations), and then the difference between them gives the true noise covariance.
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        if len(parent)!=2 or type(parent)!=type({}):
            raise Exception("science.centCovariance object should have 2 centroid parents")
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        if forGUISetup==1:
            #nsubx=self.config.getVal("wfs_nsubx")
            self.outputData=[(nsubx,nsubx,2),numpy.float64]
        else: # set up for simulation.
            #nsubx=self.config.getVal("wfs_nsubx")
            self.outputData=None#numpy.zeros((nsubx,nsubx,2),numpy.float64)

            self.sumcenta=None#numpy.zeros((nsubx,nsubx,2),numpy.float64)
            self.sumcentb=None#numpy.zeros((nsubx,nsubx,2),numpy.float64)
            self.sum2centa=None#numpy.zeros((nsubx,nsubx,2),numpy.float64)
            self.sum2centb=None#numpy.zeros((nsubx,nsubx,2),numpy.float64)
            #the covariances
            self.cova=None#numpy.zeros((nsubx,nsubx,2),numpy.float64)
            self.covb=None#numpy.zeros((nsubx,nsubx,2),numpy.float64)
            self.naverage=0
            self.control={"compute":0}
            self.lastCompute=0

      
    def newParent(self,parent,idstr=None):
        """Add a new parent."""
        if idstr not in self.idstr:
            raise Exception("xinterp_dm - newParent() idstr not known. %s %s"%(idstr,self.idstr))
        indx=self.idstr.index(idstr)
        self.parent=parent
        
    def generateNext(self,ms=None):
        if self.generate==1:
            if self.newDataWaiting:
                self.dataValid=1
                if self.control["compute"]:
                    if self.lastCompute==0:
                        self.naverage=0
                    keys=self.parent.keys()
                    keys.sort()
                    self.keya=keys[0]
                    self.keyb=keys[1]
                    if self.outputData==None:
                        self.outputData=numpy.zeros(self.parent[self.keya].outputData.shape,numpy.float32)
                        self.sumcenta=numpy.zeros(self.parent[self.keya].outputData.shape,numpy.float32)
                        self.sumcentb=numpy.zeros(self.parent[self.keya].outputData.shape,numpy.float32)
                        self.sum2centa=numpy.zeros(self.parent[self.keya].outputData.shape,numpy.float32)
                        self.sum2centb=numpy.zeros(self.parent[self.keya].outputData.shape,numpy.float32)
                        self.cova=numpy.zeros(self.parent[self.keya].outputData.shape,numpy.float32)
                        self.covb=numpy.zeros(self.parent[self.keya].outputData.shape,numpy.float32)
                    if self.naverage==0:
                        self.sumcenta[:]=0
                        self.sumcentb[:]=0
                        self.sum2centa[:]=0
                        self.sum2centb[:]=0
                    self.sumcenta+=self.parent[self.keya].outputData
                    self.sum2centa+=self.parent[self.keya].outputData*self.parent[self.keya].outputData
                    self.sumcentb+=self.parent[self.keyb].outputData
                    self.sum2centb+=self.parent[self.keyb].outputData*self.parent[self.keyb].outputData
                    self.naverage+=1
                    self.cova[:]=self.sum2centa/self.naverage-self.sumcenta*self.sumcenta/(self.naverage*self.naverage)
                    self.covb[:]=self.sum2centb/self.naverage-self.sumcentb*self.sumcentb/(self.naverage*self.naverage)
                    self.outputData[:]=self.covb-self.cova
                    if numpy.sum(self.outputData)<0:
                        self.outputData*=-1
                self.lastCompute=self.control["compute"]
        else:
            self.dataValid=0
            
    def postCovariances(self):
        """Post the covariance for other objects to use:
        they can use config.postGet("centroidNoiseCovariance")."""
        if self.idstr[0]!=None:
            name="centroidNoiseCovariance_"+self.idstr[0]
        else:
            name="centroidNoiseCovariance"
        self.config.post(name,self.outputData)
        util.FITS.Write(self.outputData,"centroidNoiseCovariance.fits")


    
    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        
        if self.idstr[0]==None or self.idstr[0]=="":
            id=""
        else:
            id=" (%s)"%self.idstr[0]
        txt=""
        txt+="""<plot title="centCov a%s" cmd="data=%s.cova" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        txt+="""<plot title="centCov b%s" cmd="data=%s.covb" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        txt+="""<plot title="centCov compute%s" when="cmd" ret="feedback" texttype="1" wintype="mainwindow">\n<cmd>\nfeedback=1-%s.control['compute']\n%s.control['compute']=feedback\nif feedback==0:\n %s.postCovariances()\n</cmd>\nbutton=feedback\nif feedback==0:\n msg='Noise covariance finished computing'\nelse:\n msg='Computing noise covariance...'\nfeedback=msg</plot>\n"""%(id,objname,objname,objname)
        txt+="""<plot title="centCov noise covariance%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt"/>\n"""%(id,objname)
        return txt


