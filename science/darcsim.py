import os
import tempfile
import socket
import traceback
import numpy
import base.aobase
import util.FITS,util.dist
import util.zernikeMod
import darc

class Darc(base.aobase.aobase):
    """
    Simulation object that takes wfs images for input, and provides dm demands as output.  darc is used for intermediate calculations...

What info do I need:
ncam, npxlx,npxly, nacts
subaplocation.  nsubaps -> subapFlag.
ncamThreads
bgimage
darknoise, flatfield (both None for now)
rmx
gain
refCentroids?
reconstructMode = truth
thresholdValue, thresholdType
E (need to modify darc to allow this to be None).
decayFactor.

Also an option to allow an existing darc to connect/disconnect/reconnect to the simulation.  To allow swapping between physical/simulation.  The swapping should be done by darc commands, outside of the simulation.  e.g. setting the mirror library and camera library to socket versions, etc.

    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        if type(parent)!=type({}):
            parent={"1":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.pupil=self.config.getVal("pupil")
        self.atmosGeom=self.config.getVal("atmosGeom")
        self.dmObj=self.config.getVal("dmObj")
        self.dmList=self.dmObj.makeDMList(self.idstr[0])
        self.useExistingDarc=self.config.getVal("useExistingDarc",default=0)
        if len(self.dmList)==0:
            raise Exception("No DMs found - do you need to specify actuatorsFrom='%s' in your config file?"%str(self.idstr))
        self.nacts=0
        self.nactsList=[]
        self.nactsCumList=[0]#cumulative version.
        self.closedLoopList=[]
        self.dmPupList=[]
        self.npokesList=[]
        self.npokesCumList=[0]
        for dm in self.dmList:
            if dm.zonalDM:
                tmp=dm.computeDMPupil(self.atmosGeom,centObscuration=self.pupil.r2,retPupil=0)
                # tmp is dmflag,subarea (or None,None for modal DMs.)
                self.dmPupList.append(tmp[0])
            
                self.nactsList.append(int(numpy.sum(tmp[0].ravel())))
                if dm.pokeSpacing!=None:
                    self.npokesList.append(dm.pokeSpacing**2)
                else:
                    self.npokesList.append(self.nactsList[-1])
                self.npokesCumList.append(self.npokesCumList[-1]+self.npokesList[-1])
            else:#a modal DM
                self.dmPupList.append(None)
                self.nactsList.append(dm.nact)#nact is the number of modes
                self.npokesList.append(dm.nact)
                self.npokesCumList.append(self.npokesCumList[-1]+dm.nact)
            self.nactsCumList.append(self.nactsList[-1]+self.nactsCumList[-1])
            self.closedLoopList.append(dm.closedLoop)
            
        self.nacts=self.config.getVal("nactsDarc",default=sum(self.nactsList))
        
        if forGUISetup==1:
            self.outputData=[(self.nacts,),numpy.float32]
        else:
            self.outputDataBuffer=numpy.zeros((self.nacts*4+8,),numpy.uint8)#8 byte header.
            self.outputData=self.outputDataBuffer[8:].view(numpy.float32)
            self.minarea=self.config.getVal("wfs_minarea")
            self.pokeActMapFifo=[]
            print "darcsim: Using %d DMs for reconstruction."%len(self.dmList)
            self.ngsList=self.atmosGeom.makeNGSList(self.idstr[0],minarea=self.minarea)#self.nsubxDict,None)
            self.lgsList=self.atmosGeom.makeLGSList(self.idstr[0],minarea=self.minarea)
            self.ncents=0
            self.ncentList=[]
            indiceList=[]
            for gs in self.ngsList+self.lgsList:
                subflag=self.pupil.getSubapFlag(gs.nsubx,gs.minarea)
                indiceList.append(numpy.nonzero(subflag.ravel())[0])
                self.ncentList.append(numpy.sum(subflag.ravel()))
                #self.ncents+=2*gs.nsubx**2
            self.nwfs=len(indiceList)
            self.ncents=sum(self.ncentList)*2
            if self.ncents==0:
                raise Exception("No centroids found for darcsim %s - check that your atmos.source objects contain a reconList=['%s'] argument"%(idstr,idstr))
            self.centIndex=numpy.zeros((self.ncents,),numpy.int32)
            pos=0
            for i in xrange(len(self.ncentList)):
                self.centIndex[pos:pos+self.ncentList[i]]=(indiceList[i]*2).astype(numpy.int32)
                self.centIndex[pos+self.ncents/2:pos+self.ncents/2+self.ncentList[i]]=(indiceList[i]*2+1).astype(numpy.int32)
                pos+=self.ncentList[i]
            self.wfsIDList=[]
            wfsIDList=self.atmosGeom.getWFSOrder(self.idstr[0])#self.config.getVal("wfsIDList",default=self.parent.keys())#the list of parent objects (specifies the order in which the centroids are placed - if that makes a difference?)...
            #print wfsIDList
            #Note:  If the parent is a single module that has already combined all the WFS information, this is fine.  It should have a wfscent_X value for npxlx of the number of pixels, and npxly of 1.  Everything else won't be used.
            for key in wfsIDList:
                if key not in self.parent.keys():
                    #if get this error, could try specifying a dummy parent that has outputData an array of zeros.
                    raise Exception("darcsim: key %s not found in parent, so not using"%str(key))
                else:
                    self.wfsIDList.append(key)
            #The length of wfsIDList is ncam.
            ncam=len(self.wfsIDList)
            #Now create npxlx and npxly lists...
            npxlx=[]
            npxly=[]
            so=self.config.searchOrder[:]
            sl=[]
            nsubList=[]
            print "Wavefront sensors for darcsim %s: %s"%(str(self.idstr),str(self.wfsIDList))
            for key in self.wfsIDList:
                self.config.setSearchOrder(["wfscent_%s"%key,"wfscent","globals"])
                pupfn=self.config.getVal("pupil")
                if type(pupfn)!=numpy.ndarray:
                    pupfn=pupfn.fn
                pupfn=pupfn.astype(numpy.float32)
                wfs_minarea=self.config.getVal("wfs_minarea")
                nsubx=self.config.getVal("wfs_nsubx")
                n=self.config.getVal("wfs_n")
                nimg=self.config.getVal("wfs_nimg")
                y,x=nsubx*nimg,nsubx*nimg
                npxlx.append(x)
                npxly.append(y)
                nsub=0
                for i in xrange(nsubx):        
                    for j in xrange(nsubx):
                        if pupfn[i*n:(i+1)*n,j*n:(j+1)*n].sum()>wfs_minarea*n*n:
                            #this is a used subap.
                            sl.append((i*nimg,(i+1)*nimg,1,j*nimg,(j+1)*nimg,1))
                            nsub+=1
                nsubList.append(nsub)
            if len(self.wfsIDList)==0:#parent is probably something that combines wfs images.
                npxlx=self.config.getVal("npxlx")
                npxly=self.config.getVal("npxly")
            self.config.setSearchOrder(so)
            #Open a listening socket, for the darc camera and mirror to connect to.
            self.lsock=socket.socket()
            self.lsock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            port=8000
            while port<9000:
                print "binding to port %d..."%port
                try:
                    self.lsock.bind(("",port))
                except:
                    port+=1
                else:
                    break
            if port==9000:
                raise Exception("Unable to bind...")
            print "Bound to port %d"%port
            self.port=port
            self.lsock.listen(2)
            self.frameno=0

            if self.useExistingDarc:
                
                npxls=(numpy.array(npxlx)*numpy.array(npxly)).sum()
                npxls=self.config.getVal("npxlsDarc",npxls)
                self.npxls=npxls
                print "Expecting to deliver %d pixels to darc"%npxls
                self.inputData=numpy.zeros((npxls+2,),numpy.float32)
                self.inputHdr=self.inputData[:2].view(numpy.uint32)
            else:#start darc for this simulation.
                fname=self.config.getVal("darcConfigFileName",default=None,raiseerror=0)
                if fname==None or not os.path.exists(fname) or self.config.getVal("overwriteDarcConfig",default=1):
                    #need to create a new config...
                    subapLocation=numpy.array(sl).astype(numpy.int32)
                    nslopes=subapLocation.shape[0]*2
                    bg=self.config.getVal("darcBG",default=0.)
                    rmx=self.config.getVal("darcRmx",default=0.)
                    gain=self.config.getVal("darcGain",default=0.2)
                    decay=self.config.getVal("darcDecay",default=0.99)
                    E=self.config.getVal("darcEmx",default=None,raiseerror=0)
                    nthreads=self.config.getVal("darcThreads",default=1)
                    thrVal=self.config.getVal("darcThresholdValue",default=0.)
                    thrType=self.config.getVal("darcThresholdType",default=1)
                    self.darcArgs=self.config.getVal("darcArgs",default="")
                    self.prefix=self.config.getVal("darcPrefix",default="sim%s"%self.idstr[0])
                    self.prefix+="b%d"%self.config.batchno
                    npxls=(numpy.array(npxlx)*numpy.array(npxly)).sum()
                    self.npxls=npxls
                    self.inputData=numpy.zeros((npxls+2,),numpy.float32)
                    self.inputHdr=self.inputData[:2].view(numpy.uint32)
                    # now make the darc config file...
                    sltxt=arrToDarcTxt(subapLocation)
                    bgtxt=self.valToDarcArrTxt(bg,npxls)
                    rmxtxt=self.valToDarcArrTxt(rmx,(self.nacts,nslopes))
                    gaintxt=self.valToDarcArrTxt(gain,self.nacts)
                    decaytxt=self.valToDarcArrTxt(decay,self.nacts)
                    nthreadstxt=self.valToDarcArrTxt(nthreads,len(self.wfsIDList))
                    darcUserTxt=self.config.getVal("darcUserTxt",default="")
                    txt="""
import numpy
cameraParams=numpy.zeros((5,),numpy.int32)
cameraParams[0]=1#asfloat
cameraParams[1]=%d#port
cameraParams[2:]=numpy.fromstring("127.0.0.1\\0\\0\\0",dtype=numpy.int32)
mirrorParams=numpy.zeros((10,),numpy.int32)
mirrorParams[0]=1#timeout
mirrorParams[1]=%d#port
mirrorParams[2]=1#affin el size
mirrorParams[3]=1#priority
mirrorParams[4]=-1#affin
mirrorParams[5]=0#sendprefix
mirrorParams[6]=1#asfloat
mirrorParams[7:]=numpy.fromstring("127.0.0.1\\0\\0\\0",dtype=numpy.int32)

control={"ncam":%d,
"nsub":%s,
"npxlx":%s,
"npxly":%s,
"nacts":%d,
"ncamThreads":%s,
"subapLocation":%s,
"bgImage":%s,
"rmx":%s,
"gain":%s,
"decayFactor":%s,
"E":None,
"reconstructMode":"truth",
"thresholdValue":%g,
"thresholdType":%d,
"cameraName":"libcamsocket.so",
"cameraParams":cameraParams,
"camerasOpen":1,
"mirrorName":"libmirrorSocket.so",
"mirrorOpen":1,
"mirrorParams":mirrorParams,
"bleedGain":%g,
"v0":numpy.zeros((%d,),numpy.float32),
"actMin":-numpy.ones((%d,),numpy.float32)*1e6,
"actMax":numpy.ones((%d,),numpy.float32)*1e6,
}
%s
"""%(self.port,self.port,ncam,str(nsubList),str(npxlx),str(npxly),self.nacts,nthreadstxt,sltxt,bgtxt,rmxtxt,gaintxt,decaytxt,thrVal,thrType,0.01/self.nacts,self.nacts,self.nacts,self.nacts,darcUserTxt)
                    if self.config.getVal("printDarcConfigFile",default=1):
                        print txt
                    if fname==None:
                        fd=tempfile.NamedTemporaryFile(mode="w",prefix="config",suffix=".py")
                        fname=fd.name
                    else:
                        fd=open(fname,"w")
                    fd.write(txt)
                    fd.close()
                    print "Written config to %s"%fname
                else:
                    print "Using darc config file %s"%fname
                self.fname=fname
            #now start darc, and wait for it to connect.
            if not self.useExistingDarc:
                self.startDarc()
            print "Waiting for darc to connect"
            self.camsock,addr=self.lsock.accept()
            print "Accepted camera connection from %s"%str(addr)
            self.mirsock,addr=self.lsock.accept()
            print "Accepted mirror connection from %s"%str(addr)
            #self.ctrl=darc.Control(self.prefix)

    def arrStrEvalable(self,arr):
        """Converts arr to a string that will return arr if evaled"""
        txt=str(arr)
        go=1
        while go==1:
            try:
                if txt.index("  "):
                    txt=txt.replace("  "," ")
            except:
                go=0
        txt="numpy.array("+txt.replace("[ ","[").replace(" ",",")+")"
        return txt
        

    def valToDarcArrTxt(self,val,size,dtype="numpy.float32"):
        if type(size) not in (type([]),type(())):
            size=(size,)
        if type(val) in [type(0),type(0.)]:
            txt="numpy.ones(%s,%s)*%g"%(size,dtype,val)
        elif type(val)==None:
            txt="None"
        elif type(val)==type(""):
            if os.path.exists(val):
                txt="""FITS.Read("%s")[1]"""%val
            else:
                txt=val
        elif type(val)==numpy.ndarray:
            txt=self.arrStrEvalable(val)
        else:
            raise Exception("Unknown value %s"%str(val))
        return txt

    def newParent(self,parent,idstr=None):
        raise Exception("Please don't call newParent for darcsim module")
        
    def generateNext(self,ms=None):
        if self.generate==1:
            if self.newDataWaiting:
                self.dataValid=1
                for key in self.wfsIDList:
                    if self.parent[key].dataValid==0:
                        self.dataValid=0
            if self.dataValid:
                self.getInputData()
                self.sendToDarc()
                self.getDarcOutput()
        else:
            self.dataValid=0
            

    def startDarc(self):
        cmd="darccontrol %s --prefix=%s %s &"%(self.fname,self.prefix,self.darcArgs)
        print "Running: %s"%cmd
        os.system(cmd)

    def reopenSockets(self):
        """Reopens the listening sockets and waits for darc to connect"""
        print "Waiting for darc to reconnect cameras and mirrors... "
        self.camsock,addr=self.lsock.accept()
        print "Accepted camera connection from %s"%str(addr)
        self.mirsock,addr=self.lsock.accept()
        print "Accepted mirror connection from %s"%str(addr)

    def __del__(self):
        self.endSim()
    def endSim(self):
        if not self.useExistingDarc:
            print "Stopping darc %s"%self.prefix
            try:
                self.ctrl=darc.Control(self.prefix)
                if self.ctrl!=None:
                    self.ctrl.ControlHalt()
                    self.ctrl=None
            except:
                print "Unable to halt darc %s"%self.prefix
                traceback.print_exc()
    def getInputData(self):
        pos=2
        for key in self.wfsIDList:
            size=self.parent[key].outputData.size
            self.inputData[pos:pos+size]=self.parent[key].outputData.ravel()
            pos+=size

    def sendToDarc(self):
        self.inputHdr[0]=0xa<<28|self.npxls
        self.inputHdr[1]=self.frameno
        try:
            self.camsock.sendall(self.inputData)
        except:
            if self.useExistingDarc:
                #darc has swapped to a different camera library - so here, just listen again.
                print "Socket error in darcsim - reopening"
                self.reopenSockets()#including the mirror socket.  
                #and send again.
                print "darcsim:  Resending data"
                self.camsock.sendall(self.inputData)
            else:#something went wrong.
                raise

    def getDarcOutput(self):
        nrec=0
        while nrec<self.outputDataBuffer.size:
            r=self.mirsock.recv_into(self.outputDataBuffer[nrec:])
            if r>0:
                nrec+=r
            else:
                if self.useExistingDarc:
                    print "Socket error receiving DM demands - reopening sockets"
                    self.reopenSockets()
                else:
                    raise Exception("Error receiving... %d (received %d)"%(r,nrec))
    
    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr==None:
            id=""
        else:
            id=" (%s)"%self.idstr
        txt=""
        txt+="""<plot title="darcsim output%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
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
        paramList.append(base.dataType.dataType(description="telDiam",typ="f",val="42.",comment="TODO: Telescope diameter (m)"))
        paramList.append(base.dataType.dataType(description="npup",typ="i",val="1024",comment="TODO: Number of pixels used to sample the pupil"))
        paramList.append(base.dataType.dataType(description="nAct",typ="eval",val="wfs_nsubx+1",comment="TODO: Number of actuators across the pupil"))
        paramList.append(base.dataType.dataType(description="wfs_n",typ="eval",val="npup/wfs_nsubx",comment="TODO: Number of pixels across a subap"))
        paramList.append(base.dataType.dataType(description="atmosGeom",typ="code",val="import util.atmos;atmosGeom=util.atmos.geom(layerDict, sourceList,ntel,npup,telDiam)",comment="TODO: atmosGeom with arguments layerDict, sourceList,ntel,npup,telDiam.  layerDict is a dictionar with keys equal to the layer name (idstr used by infScrn object) and values equal to a tuple of (height, direction, speed, strength, initSeed), and sourceList is a list of ources equal to a tuple of (idstr (infAtmos), theta, phi, alt, nsubx or None)."))
        paramList.append(base.dataType.dataType(description="dmObj",typ="code",val="import util.dm;dmObj=util.dm.dmOverview(dmInfoList,atmosGeom=atmosGeom)",comment="TODO: dmInfoList is a list of util.dm.dmInfo objects, initialised with (label (for this dm), idlist (list of (dm ID,source ID) or just sourceID, the idstr for a particular DM object (height and direction) and the idstr for a given source direction), height, nact, and various other things for which the defaults may be okay (see the source code))"))
        
        paramList.append(base.dataType.dataType(description="pupil",typ="code",val="import util.tel;pupil=util.tel.Pupil(npup,npup/2,npup/2*telSec/telDiam,wfs_nsubx,wfs_minarea)",comment="Telescope pupil object"))
        return paramList
	
def arrToDarcTxt(arr):
    numpy.set_printoptions(threshold=arr.size+1)
    txt=str(arr)
    numpy.set_printoptions(threshold=1000)
    go=1
    while go==1:
        try:
            if txt.index("  "):
                txt=txt.replace("  "," ")
        except:
            go=0
    txt="numpy.array("+txt.replace("[ ","[").replace(" ",",")+")"
    return txt
