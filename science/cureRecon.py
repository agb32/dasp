import numpy
import ccure
import base.aobase
import tel

class recon(base.aobase.aobase):

    def __init__(self, parent, config, args={}, forGUISetup=0, debug=None, idstr=None):
        print("**** CureRecon in use ****")
        global util
        if type(parent)!=type({}):
            parent={"cent":parent}
        base.aobase.aobase.__init__(self, parent, config, args, forGUISetup=forGUISetup, debug=debug, idstr=idstr)
        self.dataValid=1#data is valid even before the first iteration because we assume the mirror starts zerod
        self.pupil=self.config.getVal("pupil")
        self.cureScale=self.config.getVal("cureScale", default=1.0)
        self.atmosGeom=self.config.getVal("atmosGeom")
        self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
        if self.dmObj==None or type(self.dmObj)!=type(self.atmosGeom):
            print "DEPRECATION: warning: dmObj should now be dmOverview"
            self.dmObj=self.config.getVal("dmObj")
        self.dmList=self.dmObj.makeDMList(self.idstr[0])
        if len(self.dmList)==0:
            raise Exception("No DMs found - do you need to specify actuatorsFrom='%s' in your config file?"%str(
                    self.idstr))

        self.nacts=0
        self.nactsList=[]
        self.nactsCumList=[0]#cumulative version.
        self.closedLoopList=[]
        self.dmPupList=[]
        self.npokesList=[]
        self.npokesCumList=[0]
        self.reconDtype=self.config.getVal("reconDtype", default=numpy.float32)
        self.dmModeType=self.config.getVal("dmModeType", default="poke")

        self.nsub=self.config.getVal("wfs_nsubx")
        self.minarea=self.config.getVal("wfs_minarea",raiseerror=0)
        self.ngsList=self.atmosGeom.makeNGSList(self.idstr[0], minarea=self.minarea, pupil=self.pupil)
        self.lgsList=self.atmosGeom.makeLGSList(self.idstr[0], minarea=self.minarea, pupil=self.pupil)

        for dm in self.dmList:
            if dm.zonalDM:
                self.dmFlag=dm.getDMFlag(self.atmosGeom, centObscuration=self.pupil.r2)

                self.dmPupList.append(self.dmFlag)
                self.nactsList.append(int(self.dmFlag.sum()))
                if self.dmModeType=="poke":
                    if dm.pokeSpacing!=None:
                        self.npokesList.append(dm.pokeSpacing**2)
                    else:
                        self.npokesList.append(self.nactsList[-1])
                    self.npokesCumList.append(self.npokesCumList[-1] + self.npokesList[-1])
                elif self.dmModeType=="modalPoke":
                    self.npokesList.append(0)
            else:#a modal DM
                raise Exception("Modal DMs are not yet implemented with CuRe reconstruction")

            self.nactsCumList.append(self.nactsList[-1] + self.nactsCumList[-1])
            self.closedLoopList.append(dm.closedLoop)

        self.nacts=sum(self.nactsList)

        self.ncents=0
        self.ncentList=[]
        indiceList=[]
        for gs in self.ngsList+self.lgsList:
            subflag=gs.getSubapFlag()
            indiceList.append(numpy.nonzero(subflag.ravel())[0])
            self.ncentList.append(numpy.sum(subflag.ravel()))
        self.nwfs=len(indiceList)
        self.ncents=sum(self.ncentList)*2
        if self.ncents==0:
            raise Exception("No centroids found for CuReReon %s - check that your atmos.source objects contain a reconList=['%s'] argument"%(idstr,idstr))
        self.centIndex=numpy.zeros((self.ncents,),numpy.int32)

        pos=0
        for i in xrange(len(self.ncentList)):
            self.centIndex[pos:pos+self.ncentList[i]] = (indiceList[i]*2).astype(numpy.int32)
            self.centIndex[pos+self.ncents/2:pos + self.ncents/2 + self.ncentList[i]]=\
                (indiceList[i]*2+1).astype(numpy.int32)
            pos+=self.ncentList[i]


        flag = self.ngsList[0].getSubapFlag()
        self.phaseFlag = flag
        self.numPhasePoints = flag.sum()
        
        # Creates a square WFS flag for initialising CuRe
        squareFlag = numpy.ones(flag.shape).astype(numpy.int32)
        self.numPhasePoints = squareFlag.sum()

        # Creates a square DM flag for initialising CuRe
        self.fakeDmFlag = numpy.ones(self.dmFlag.shape).astype(numpy.int32)
        self.nacts = self.fakeDmFlag.sum()

        self.cure = ccure.open(self.nsub,self.nacts, squareFlag.astype(numpy.int32),1)
        
        self.cgain=numpy.ones((self.nacts,),"f") * self.cureScale
        self.mapObj=numpy.identity(self.nacts).astype("f")

        self.outputData=numpy.zeros((self.dmFlag.sum(),), numpy.float32)
        self.inputData=numpy.zeros(self.numPhasePoints*2, numpy.float32)

        self.wfsIDList=[]
        wfsIDList=self.atmosGeom.getWFSOrder(self.idstr[0])
        for key in wfsIDList:
            if key not in self.parent.keys():
                #if get this error, could try specifying a dummy parent that has outputData an array of zeros.
                raise Exception("CuReRecon: key %s not found in parent, so not using"%str(key))
            else:
                self.wfsIDList.append(key)
        self.parentDataValid=numpy.zeros((len(self.wfsIDList)), numpy.int32)

    def generateNext(self):
        if self.generate==1:
            if self.newDataWaiting:#this is always 1(!)
                nin=0
                for key in self.parent.keys():
                    if self.parent[key].dataValid==1:
                        nin+=1

                if nin>0:
                    if nin==len(self.parent.keys()):
                        self.dataValid=1
                    elif self.multirate==1:
                        self.dataValid=1#some WFSs are valid
                    else:
                        print("WARNING:**CuReD Recon**: got some data but not "+
                           "all, setting dataValid=0")
                        self.dataValid=0
                else:
                    self.dataValid=0
            if self.dataValid:
                self.getInput() # Parses the data into the correct form
                self.calcData()# Calculate DM updates
        else:
            self.dataValid=0

    def getInput(self):
        cnt=0
        for i in range(len(self.wfsIDList)):
            key=self.wfsIDList[i]
            ns=self.ncentList[i]

            self.parentDataValid[i]=self.parent[key].dataValid

            if self.parent[key].dataValid==1:
                self.inputData[:self.numPhasePoints] = self.parent[key].outputData[:,:,0].ravel()
                self.inputData[self.numPhasePoints:] = self.parent[key].outputData[:,:,1].ravel()
            cnt+=ns

    def calcData(self):
        # NaN values can be produced by the WFS
        # This catches them and sets them to zero
        nones = numpy.count_nonzero(numpy.isnan(self.inputData))
        if (nones > 0):
            self.inputData[numpy.flatnonzero(numpy.isnan(self.inputData))] = 0
        
        xs = self.inputData[:self.numPhasePoints].astype('f')
        ys = self.inputData[self.numPhasePoints:].astype('f')

        if (self.exportData):
            self.xList.append(xs)
            self.yList.append(ys)

        out=numpy.zeros((self.nacts,),numpy.float32)
        ccure.run(self.cure,xs,ys,out,self.cgain,self.mapObj,1)

        indexs = numpy.flatnonzero(self.dmFlag)
        self.outputData = out[indexs]

        return

    ### Called at the end of the simulation to close the CuRe module
    def endSim(self):
        ccure.close(self.cure)