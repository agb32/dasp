#$Id: el_dm.py,v 1.41 2009/06/05 12:53:11 ali Exp $

### Segmented DM simulation
### Add DM phases to input phase array and forward to output
### Update mirror figure using control signals from el_recon


import base.aobase
import math,thread
#import Numeric
import time
import numpy
import cPickle,types
from sys import argv
from util.tel import Pupil

    

class dm(base.aobase.aobase):
    """Parent should be a dictionary of recon and atmos.
    This object is able to implement resource sharing, but only to the extent
    of having one reconstructor parent - ie it is able to physically represent
    one DM, taking atmos from lots of sky directions.  
    The mirror shape itself is only updated for the first resource sharing
    object, and then stays like this for all the others.
    Note that this DM object can only have 1 unique conjugate height -
    otherwise it would be unphysical...
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        if type(parent)!=type({}):
            parent={"1":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.atmosPhaseType=self.config.getVal("atmosPhaseType",default="phaseonly")
        if self.atmosPhaseType not in ["phaseonly","phaseamp","realimag"]:
            raise Exception("el_dm: atmosPhaseType not known %s"%self.atmosPhaseType)
        if forGUISetup==1:
            npup=self.config.getVal("npup")
            if self.atmosPhaseType=="phaseonly":
                self.outputData=[(npup,npup),numpy.float64]
            elif self.atmosPhaseType=="phaseamp":
                self.outputData=[(2,npup,npup),numpy.float64]
            elif self.atmosPhaseType=="realimag":
                self.outputData=[(npup,npup),numpy.complex64]
        else:
            #self.control[zerodm and closedm] are controlled by pyro.
            self.control={"dm_update":1,"useTilt":1}#,"zero_dm":0,"close_dm":1}
            #self.phaseIsNone=0
            npup=self.npup=self.config.getVal("npup")                    # Pupil phase array size - ie number of pixels from phasescreen, and passed to output.
            nsubx=self.wfs_nsubx=self.config.getVal("wfs_nsubx")          # SH array size - actually, if using a non-ground conjugated DM, nsubx may be larger than that used for the wavefront sensors - in this case, the outputs from several centroiders will be combined by a reconstructor in the correct way to be able to shape the larger DM...
            n=self.wfs_n=self.config.getVal("wfs_n")                  # Sub-aperture phase array size
            #we now need to work out which part of the DM to use... this depends on conjugate height, and source location.
            self.pupil=self.config.getVal("pupil")
            self.atmosGeom=self.config.getVal("atmosGeom",default=None,raiseerror=0)
            self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
            if self.dmObj==None or type(self.dmObj)!=type(self.atmosGeom):
                print "DEPRECATION: warning: dmObj should now be dmOverview"
                self.dmObj=self.config.getVal("dmObj",default=None,raiseerror=0)
            if self.dmObj!=None:
                self.conjHeight=self.dmObj.getHeight(self.idstr[0])#self.config.getVal("dmConjugateHeight",default=0.)#the height at which this pupil phase should be conjugated too.  One unique height for this (resource sharing) DM object.
                self.dmpup=self.dmObj.calcdmpup(self.idstr[0])#self.config.getVal("dmpup",default=self.npup)#number of pixels to use to store the DM phase - this may be larger than npup if the DM is conjugated somewhere other than ground.
                self.subflag=self.dmObj.computeDMPupil(self.idstr[0],centObscuration=self.pupil.r2,retPupil=0)[0]
                sourceID=self.dmObj.getSourceID(idstr[0])
                if self.atmosGeom!=None:
                    self.sourceLam=self.atmosGeom.sourceLambda(sourceID)
                else:
                    self.sourceLam=None
                if self.sourceLam==None:# Wavelength for the optical path for which the mirror is shaped.
                    self.sourceLam=self.config.getVal("sourceLam")
            else:
                print "el_dm: dmObj not found in param file - continuing in depreciated mode."
                self.conjHeight=self.config.getVal("dmConjugateHeight",default=0.)#depreciated
                self.dmpup=self.config.getVal("dmpup",default=self.npup)#depreciated
                self.subflag=self.pupil.subflag
                sourceID=self.config.getVal("sourceID",default=self.idstr[0])
                self.sourceLam=self.dmObj.getDM(self.idstr[0]).reconLam
            self.xtiltfn=(numpy.fromfunction(self.tilt,(n,n))-float(n)/2.+0.5)/float(n)
            self.ytiltfn = numpy.transpose(self.xtiltfn)          # Subap tilt fns in x & y
            self.zerotilt=numpy.zeros((nsubx,nsubx),numpy.float64)
            self.xtilt=self.ytilt=self.pist=self.zerotilt
            self.phssub=numpy.zeros((n,n),numpy.float64)
            #self.pupsub=self.pupil.pupsub
            #self.subarea=self.pupil.subarea
            self.telDiam=self.config.getVal("telDiam")

          
            #wfsLam=self.config.getVal("wfslam",default=sourceLam)
            #self.wavelengthRatio=wfsLam/sourceLam#wfs_lam/lam
            
            self.dmphs=numpy.zeros((self.dmpup,self.dmpup),numpy.float64)   # DM figure (optical phase) 
            if self.atmosPhaseType=="phaseonly":
                self.outputData=numpy.zeros((npup,npup),numpy.float64)     # Input phase array
            elif self.atmosPhaseType=="phaseamp":
                self.outputData=numpy.zeros((2,npup,npup),numpy.float64)     # Input phase array
            elif self.atmosPhaseType=="realimag":
                self.outputData=numpy.zeros((npup,npup),numpy.complex64)     # Input phase array
            for i in xrange(len(self.idstr)):
                self.initialise(self.parentList[i],self.idstr[i])

    def initialise(self,parentDict,idstr):
        """note, parent may be a dictionary here (atmos and recon)..."""
        if type(parentDict)!=type({}):
            parentDict={"1":parentDict}
        this=base.aobase.resourceSharer(parentDict,self.config,idstr,self.moduleName)
        self.thisObjList.append(this)

        #for a source with theta and a conjugate height c, the
        #separation at this height will be c tan(theta).  As a
        #fraction of the phasescreen width, this will be c
        #tan(theta)/telDiam.  So, in pixels, this will be c
        #tan(theta)/telDiam * npup.  Multiply this by
        #cos/sin(phi), and you then have the correct offsets.
        this.xoff=(self.dmpup-self.npup)/2#central in the case of
        this.yoff=(self.dmpup-self.npup)/2#ground conjugated.
        if self.dmObj==None:
            this.sourceID=this.config.getVal("sourceID",default=idstr)#sourceID is the idstr for the source direction as provided for the atmosGeom object. depreciated
        else:
            this.sourceID=self.dmObj.getSourceID(idstr)

        if self.conjHeight!=0:
            print "todo: check that signs are right - and that xoff and yoff are the right way round... (el_dm - conjugated height)"
            if self.atmosGeom==None:
                this.sourceTheta=this.config.getVal("sourceTheta")/3600./180*numpy.pi#radians depreciated
                this.sourceAlt=this.config.getVal("sourceAlt")#depreciated
                this.sourcePhi=this.config.getVal("sourcePhi")*numpy.pi/180#radians depreciated
            else:
                this.sourceTheta=self.atmosGeom.sourceTheta(this.sourceID)*numpy.pi/180/3600.
                this.sourceAlt=self.atmosGeom.sourceAlt(this.sourceID)
                this.sourcePhi=self.atmosGeom.sourcePhi(this.sourceID)*numpy.pi/180
            if this.sourceAlt>0:
                raise Exception("el_dm implemetation error - not yet interpolating for non-infinite sources")
            r=self.conjHeight*numpy.tan(this.sourceTheta)/self.telDiam*self.npup
            this.xoff+=int(r*numpy.cos(this.sourcePhi)+0.5)#round
            this.yoff+=int(r*numpy.sin(this.sourcePhi)+0.5)#correctly
            #but don't bother with any interpolation here...
        # select the correct part of the dm phs.
        if this.xoff<0 or this.yoff<0 or this.xoff+self.npup>self.dmpup or this.yoff+self.npup>self.dmpup:
            raise Exception("DM pupil not large enough to hold all the conjugated sources... %d %d %d %d"%(this.xoff,this.yoff,self.npup,self.dmpup))
        this.selectedDmPhs=self.dmphs[this.yoff:this.yoff+self.npup,this.xoff:this.xoff+self.npup]
        if self.atmosGeom!=None:
            this.sourceLam=self.atmosGeom.sourceLambda(this.sourceID)
        else:
            this.sourceLam=None
        if this.sourceLam==None:
            this.sourceLam=this.config.getVal("sourceLam")          # Wavelength for this optical path
        #this.wfsLam=this.config.getVal("wfslam",default=this.sourceLam)
        #this.wavelengthRatio=this.wfsLam/this.sourceLam#wfs_lam/lam
        this.wavelengthAdjustor=self.sourceLam/this.sourceLam#this.wavelengthRatio/self.wavelengthRatio#the mirror will be shaped as for self.wavelengthRatio... if this.sourceLam is longer than self.sourceLam, radians of phase P-V will be smaller, so less change needed in wavelengthAdjustor.

    def finalInitialisation(self):
        """Just check that there aren't 2 objects with same idstr..."""
        tmp=[]
        for id in self.idstr:
            if id in tmp:
                raise Exception("el_dm - resource sharing cannot have 2 objects with same idstr")
            tmp.append(id)

    def newParent(self,parent,idstr=None):
        """Most likely - if created by simsetup GUI, the object will
        originally have None as parent, and newParent will be called
        with a dict of atmos and recon.  The question is how to handle
        this for the various children...
        Thats why the idstr can be passed.  The idstr is used to find the correct resource sharing object, and
        then it is given the parents.  Note, that this module will not allow 2 resource sharing objects to have
        the same idstr.
        """
        if idstr not in self.idstr:
            raise Exception("el_dm - newParent() idstr not known. %s %s"%(idstr,self.idstr))
        indx=self.idstr.index(idstr)
        self.thisObjList[indx].parent=parent

        

##             wfs_minarea=self.config.getVal("wfs_minarea")
##             for i in range(nsubx):        
##                 for j in range(nsubx):
##                     self.pupsub[i][j]=self.pupil.fn[i*n:(i+1)*n,j*n:(j+1)*n]    # Get pupil fn over subaps
##                     self.subarea[i][j]=numpy.sum(numpy.sum(self.pupsub[i][j]))                          
##                     if(self.subarea[i][j]>(wfs_minarea*n*n)):       # Flag vignetted subaps 
##                         self.subflag[i][j]=1

##             self.dmpupil=numpy.zeros((npup,npup),numpy.Float64) # DM pupil - set to zero for 
##             for i in range(nsubx):                                  # flagged subaps: assume they are    
##                 for j in range(nsubx):                              # tilted out of science fov
##                     if(self.subflag[i][j]==1):
##                         self.dmpupil[i*n:(i+1)*n,j*n:(j+1)*n]=1.
        # print thread.get_ident(),"TODO:  Should calculate reflected phase happen b4 or after dm update?"
        
##     def needData(self,msg=None):
##         return True
        
##     def getData(self,msg=None):                                 # Parent is a splitter connected to atmos and recon
##         """Parent is a splitter connected to atmos and recon"""

        
##         self.inputData=self.parent.next(self.objID,msg=msg)
##         if type(self.inputData["atmos"])==types.NoneType:
##             self.phaseIsNone=1
##         else:
##             self.phaseIsNone=0
##             self.outputData[:,]=self.inputData["atmos"]                 # Make a copy because we change phs.

            
##         reconData=self.inputData["recon"]                       # DM control data from the reconstructor 
##         if type(reconData)!=types.NoneType:                     # There was something in the feedback, or...
##             self.reconData=reconData
##             self.control["dm_update"]=1
##         else:                                                   # no feedback yet...
##             self.control["dm_update"]=0

            
    def generateNext(self,msg=None):                            # DM main loop function
        """DM main loop function"""
        t1=time.time()
        if self.debug!=None:
            print "el_dm: Doing generateNext (debug=%s)"%str(self.debug)
        this=self.thisObjList[self.currentIdObjCnt]
        if not this.parent.has_key("atmos") and not this.parent.has_key("recon"):
            print "el_dm object assigning parents automatically"
            keylist=this.parent.keys()
            for key in keylist:
                if this.parent[key].outputData.shape==(self.npup,self.npup):#the right shape for an atmos object
                    print "el_dm parent object %s becoming atmos"%str(key)
                    this.parent["atmos"]=this.parent[key]
                else:
                    print "el_dm parent object %s becoming recon"%str(key)                    
                    this.parent["recon"]=this.parent[key]
            for key in keylist:
                del(this.parent[key])
        if not this.parent.has_key("atmos"):#maybe we're just using this for a poke matrix?
            print "xinterp_dm object no parent atmos object found.  Assuming unperturbed phase in."
            class dummyAtmosClass:
                dataValid=1
                outputData=0
            this.parent["atmos"]=dummyAtmosClass()
        
        if self.generate==1:
            if self.newDataWaiting:
                if this.parent["atmos"].dataValid==1:
                    self.outputData[:,]=this.parent["atmos"].outputData#make a copy
                    self.dataValid=1
                else:
                    print "Dm: Waiting for data from atmos, but not valid"
                    self.dataValid=0
                if self.control["dm_update"]==1 and self.currentIdObjCnt==0 and this.parent["recon"].dataValid==1:
                    #we are the first of the resource sharing objects, and have reconstructor commands...
                    self.reconData=this.parent["recon"].outputData
                    #self.control["dm_update"]=1
                    l=len(self.reconData.shape)
                    if l==3:
                        self.pist=self.reconData[0]
                        if self.control["useTilt"]:
                            self.xtilt=self.reconData[1]                        # Piston+tilt update from Recon
                            self.ytilt=self.reconData[2]
                        else:
                            self.xtilt=self.ytilt=self.zerotilt
                    elif l==2:#maybe the reconstructor only returns the pistons?
                        self.xtilt=self.ytilt=self.zerotilt
                        self.pist=self.reconData
                    elif l==1:
                        self.pist=numpy.reshape(self.reconData,(self.wfs_nsubx,self.wfs_nsubx))
                    self.update() #create self.dmphs - Update the DM figure
                    self.dataValid=1#update the output...
                    
            if self.dataValid:#atmos data is valid (and/or mirror shape has changed)
                #self.source=0                                           # Sort out source !!!!!!!!!!!!!!!!!
                #lam=self.sourceLam
                #outSource=0                                             # self.sourceDict[reconproc]
                #wfs_lam should probably be a different variable in the config file... todo
                #wfs_lam=self.sourceLam

                self.selectedDmPhs=this.selectedDmPhs
                if self.atmosPhaseType=="phaseonly":
                    if this.wavelengthAdjustor==1:
                        self.outputData+=self.selectedDmPhs            # Calculate and return the reflected phase
                    else:
                        self.outputData+=self.selectedDmPhs*this.wavelengthAdjustor
                elif self.atmosPhaseType=="phaseamp":
                    if this.wavelengthAdjustor==1:
                        self.outputData+=self.selectedDmPhs            # Calculate and return the reflected phase
                    else:
                        self.outputData+=self.selectedDmPhs*this.wavelengthAdjustor
                elif self.atmosPhaseType=="realimag":
                    if this.wavelengthAdjustor==1:
                        phase=numpy.arctan2(self.outputData.imag,self.outputData.real)+self.selectedDmPhs            # Calculate reflected phase
                    else:
                        phase=numpy.arctan2(self.outputData.imag,self.outputData.real)+self.selectedDmPhs*this.wavelengthAdjustor            # Calculate reflected phase
                    amp=numpy.absolute(self.outputData)
                    self.outputData[:,]=amp*numpy.cos(phase)+1j*amp*numpy.sin(phase)
        else:
            self.dataValid=0
        if self.debug!=None:
            print "el_dm: Done generateNext (debug=%s)"%str(self.debug)
        self.generateNextTime=time.time()-t1

    def update(self):                                           # Add new pistons and tilts to mirror
        """Add new pistons and tilts to mirror"""
        nsubx=    self.wfs_nsubx
        wfs_n=    self.wfs_n
        phssub=self.phssub
##        gamma=self.gamma
##         if self.control['zero_dm']:
##             print "Zeroing DM"
##             self.control['zero_dm']=0
##             self.zxtilt*=0.0                                    # Zero the segment pistons and tilts
##             self.zytilt*=0.0
##             self.zpist*=0.0
 
        #if self.control['close_dm']:#update the DM - this is now decided in el_recon...
        for i in range(nsubx):
            for j in range(nsubx):
                if(self.subflag[i,j]==1):
                    #self.zxtilt[i][j] = self.zxtilt[i][j] + gamma*self.xtilt[i][j]# Update segment piston + tilts
                    #self.zytilt[i][j] = self.zytilt[i][j] + gamma*self.ytilt[i][j]
                    #self.zpist[i][j] =  self.zpist[i][j]  + gamma*self.pist[i][j]


                    phssub[:,]=self.xtilt[i,j]*self.xtiltfn   # Get segment phase fn for updated 
                    phssub+=self.ytilt[i,j]*self.ytiltfn  # piston and tilt
                    phssub+=self.pist[i,j]
                    self.dmphs[i*wfs_n:(i+1)*wfs_n,j*wfs_n:(j+1)*wfs_n]=phssub#*self.wavelengthRatio # Tessalate onto overall DM phs map



    def tilt(self,x,y):
	"""Create a tilt fn (z=y)"""
	return y

    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        this=self.thisObjList[self.sentPlotsCnt]
        if this.idstr==None or this.idstr=="":
            id=""
        else:
            id=" (%s)"%this.idstr
        txt=""
        if self.sentPlotsCnt==0:
            txt+="""<plot title="el_DM output%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="el_DM mirror%s" cmd="data=-%s.dmphs" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
            txt+="""<plot title="el_DM SOR pistons%s" cmd="data=-%s.pist" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        txt+="""<plot title="el_DM selected phs%s" cmd="data=-%s.thisObjList[%d].selectedDmPhs" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname,self.sentPlotsCnt)
        self.sentPlotsCnt=(self.sentPlotsCnt+1)%len(self.thisObjList)
        return txt

    def getParams(self):
        """parameters required for this module, in the form of {"paramName":defaultValue,...}
        These params can then be placed in the config file... if not set by the
        user, the param should still be in config file as default value for
        future reference purposes.
        """
        #This is a working example.  Please feel free to change the parameters
        #required. (if you do, also change the config.getParam() calls too).
        paramList=[]
        paramList.append(base.dataType.dataType(description="atmosPhaseType",typ="s",val="phaseonly",comment="Probably phaseonly"))
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="wfs_nsubx",typ="i",val="32",comment="TODO: Number of subaps across the pupil"))
        paramList.append(base.dataType.dataType(description="wfs_n",typ="eval",val="npup/wfs_nsubx",comment="TODO: Number of pixels across a subap"))
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))
        paramList.append(base.dataType.dataType(description="dmObj",typ="code",val="import util.dm;dmObj=util.dm.dmOverview(dmInfoList,atmosGeom=atmosGeom)",comment="TODO: dmInfoList is a list of util.dm.dmInfo objects, initialised with (label (for this dm), idlist (list of (dm ID,source ID) or just sourceID, the idstr for a particular DM object (height and direction) and the idstr for a given source direction), height, nact, and various other things for which the defaults may be okay (see the source code))"))
        
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(npup,npup/2,npup/2*telSec/telDiam,wfs_nsubx",comment="Telescope pupil object"))
        #paramList.append(base.dataType.dataType(description="sourceLam",typ="f",val="1650.",comment="source wavelength"))
        #paramList.append(base.dataType.dataType(description="wfslam",typ="f",val="1650.",comment="source wavelength"))

        return paramList
	

if __name__=="__main__":
    raise Exception("el_dm - cannot run standalone...")
    #warning - won't work until parentList has instances of [atmos, recon]
