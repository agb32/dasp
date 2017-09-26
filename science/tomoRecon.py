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
import os.path
import numpy
import base.aobase
import util.centroid
import util.phaseCovariance
import util.createPokeMx
import util.tomofdpcg
import util.FITS
import util.blockMatrix
import util.zernikeMod
import util.regularisation
#import cmod.svd #removed from dasp because probably depreciated.
import cmod.utils
#import util.dot as quick
quick=numpy

try:
    import scipy.linsolve
except:
    print("INFORMATION:**tomoRecon**:TO DO: Sort out import of scipy.linsolve in tomoRecon")
import scipy.sparse,scipy.linalg
import util.spmatrix
import time,types
#import Scientific.MPI
class recon(base.aobase.aobase):
    """A reconstructor for tomographic reconstruction.
    Sends actuator values for model DMs to the children.

    Currently can use the fdpcg reconstructor (though this doesn't
    close the loop!), and a sparse matrix based reconstructor (zonal
    poke matrix, MVM type of thing), a svd reconstructor, and a generalized
    inverse system (which uses SVD).

    For the spmx reconstructor, we are trying to solve Ax=b where A is the poke matrix, b is the centroids,
    and x is the reconstructed phase.  So, we first multiply by A^T (A transposed), to give A^TAx=A^Tb.  The
    matrix A^TA is then square, and so a solution for x can be found using LU decomposition.  
    Note this didn't work (not modal).  So instead, a sparse SVD implementation has been created, which should be
    able to handle the EAGLE case.
    """
    
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        global util
        if type(parent)!=type({}):
            parent={"cent":parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.dataValid=1#data is valid even before the first iteration because we assume the mirror starts zerod
        self.pupil=self.config.getVal("pupil")
        self.atmosGeom=self.config.getVal("atmosGeom")
        self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
        if self.dmObj==None or type(self.dmObj)!=type(self.atmosGeom):
            print "DEPRECATION: warning: dmObj should now be dmOverview"
            self.dmObj=self.config.getVal("dmObj")
        self.dmList=self.dmObj.makeDMList(self.idstr[0])
        self.reconObj=None#to aid addition of generic algorithms.  This must implement a few methods.  In 2014 only used by fewha, to be extended in the future
        if len(self.dmList)==0:
            raise Exception("No DMs found - do you need to specify actuatorsFrom='%s' in your config file?"%str(
                    self.idstr))
        #self.nactDMList=self.config.getVal("nactDMList")#list of nact for the DMs.
        #self.DMConjHeightList=self.config.getVal("DMConjHeightList")#list of conjugate heights for the DMs.
        #self.dmList=self.atmosGeom.makeDMList(self.nactDMList,self.DMConjHeightList)
        #self.nsubxDict=self.config.getVal("nsubxDict")#dict of idstr,nsubx for the wavefront sensors.
        self.nacts=0
        self.nactsList=[]
        self.nactsCumList=[0]#cumulative version.
        self.closedLoopList=[]
        self.dmPupList=[]
        self.npokesList=[]
        self.npokesCumList=[0]
        self.reconDtype=self.config.getVal("reconDtype",default=numpy.float32)
        self.dmModeType=self.config.getVal("dmModeType",default="poke")# specifies what the modes of
        # the dm are - are they just pokes, or a more complicated shape fitted to the actuators?
        # Note, that zernike DMs can also be poked... at the moment, only poke is supported.
        if self.dmModeType not in ["poke","modalPoke"]:
            raise Exception("tomoRecon: dmModeType should be one of: poke, modalPoke")

        for dm in self.dmList:
            if dm.zonalDM:
                tmp=dm.getDMFlag(self.atmosGeom,centObscuration=self.pupil.r2)
                self.dmPupList.append(tmp)
                self.nactsList.append(int(tmp.sum()))
                if self.dmModeType=="poke":
                    if dm.pokeSpacing!=None:
                        self.npokesList.append(dm.pokeSpacing**2)
                    else:
                        self.npokesList.append(self.nactsList[-1])
                    self.npokesCumList.append(self.npokesCumList[-1]+self.npokesList[-1])
                elif self.dmModeType=="modalPoke":
                    self.npokesList.append(0)
            else:#a modal DM
                self.dmPupList.append(None)
                self.nactsList.append(dm.nact)#nact is the number of modes
                self.npokesList.append(dm.nact)
                self.npokesCumList.append(self.npokesCumList[-1]+dm.nact)
            self.nactsCumList.append(self.nactsList[-1]+self.nactsCumList[-1])
            self.closedLoopList.append(dm.closedLoop)

        self.nacts=sum(self.nactsList)
        self.extraActs=self.config.getVal("extraActs",default=0)# won't be sent to dm, but can be useful if 
                                                                # need to do anything clever with reconstruction.
        self.extraActsDecay=self.config.getVal("extraActsDecay",default=0)
        self.nacts+=self.extraActs
        self.npokes=sum(self.npokesList)
        if forGUISetup==1:
            self.outputData=[(self.nacts,),numpy.float64]
        else:
            self.mirrorScale=None
            self.scaleWhenPoking=self.config.getVal("scaleWhenPoking",default=0)
            self.minarea=self.config.getVal("wfs_minarea",raiseerror=0)
            self.pokeActMapFifo=[]
            print(("INFORMATION:**tomoRecon**: Using {0:d} DMs for "+
                  "reconstruction").format( len(self.dmList) ) )
##(old)            print("INFORMATION:**tomoRecon**: Using %d DMs for reconstruction."%len(self.dmList))
            self.ngsList=self.atmosGeom.makeNGSList(self.idstr[0],minarea=self.minarea,pupil=self.pupil)#self.nsubxDict,None)
            self.lgsList=self.atmosGeom.makeLGSList(self.idstr[0],minarea=self.minarea,pupil=self.pupil)
            self.subtractTipTilt=self.config.getVal("subtractTipTilt",default=0)
            self.ncents=0
            self.ncentList=[]
            indiceList=[]
            for gs in self.ngsList+self.lgsList:
                subflag=gs.getSubapFlag()#self.pupil.getSubapFlag(gs.nsubx,gs.minarea)
                indiceList.append(numpy.nonzero(subflag.ravel())[0])
                self.ncentList.append(numpy.sum(subflag.ravel()))
                #self.ncents+=2*gs.nsubx**2
            self.nwfs=len(indiceList)
            self.ncents=sum(self.ncentList)*2
            if self.ncents==0:
                raise Exception("No centroids found for tomoRecon %s - check that your atmos.source objects contain a reconList=['%s'] argument"%(idstr,idstr))
            self.centIndex=numpy.zeros((self.ncents,),numpy.int32)
            pos=0
            for i in xrange(len(self.ncentList)):
                self.centIndex[pos:pos+self.ncentList[i]]=(indiceList[i]*2).astype(numpy.int32)
                self.centIndex[pos+self.ncents/2:pos+self.ncents/2+self.ncentList[i]]=\
                    (indiceList[i]*2+1).astype(numpy.int32)
                pos+=self.ncentList[i]
            self.reconType=self.config.getVal("recontype",default="pcg")
            supportedReconTypes = ["pcg","fdpcg","spmx","spmxSVD","spmxGI","svd","SVD",
                                   "MAP","pinv","reg","regSmall","regBig","regularised","dicure","fewha"]
            #not sure what this does.
            #if "newReconType" in dir():
            #    supportedReconTypes.append(newReconType)
            # check if the recontype given in the parameter file is valid:
            if self.reconType not in supportedReconTypes:
                raise ValueError("tomoRecon: recontype must be one of:", supportedReconTypes)
            del supportedReconTypes

            if self.reconType=="spmx":
                raise Exception("tomoRecon - reconType spmx does not work!")
            if self.reconType in ["spmx","spmxSVD","spmxGI"]:
                self.spmx=None
                self.sparsePmxType=self.config.getVal("sparsePmxType",default="lil")#whether to 
                # use scipy.sparse.lil_matrix or svd.sparseMatrixCreate... or a mmap'd dense
                # matrix which is then sparsified.  Can be one of "lil", "csc","mmap" with mmap
                # being fastest, though requiring space fora mmap file.
                if self.sparsePmxType not in ["lil","csc","dense"]:
                    print("WARNING:**tomoRecon**:: sparsePmxType not known - "+
                        "will assume dense but will not sparsify")
                if self.sparsePmxType=="csc":
                    # default 10% sparsity:
                    self.ndata=int(self.config.getVal("spmxNdata",default=self.nmodes*self.ncents*0.1))
                    if self.ndata>=1<<31:
                        self.ndata=(1<<31)-1
                        print(("ERROR:**tomoRecon**: spmxNdata>=1<<31 - must "+
                              "fit into an int... downsizing to "+
                              str(self.ndata)))

            if self.reconType=="SVD":
                self.reconType="svd"
            self.npup=config.getVal("npup")
            self.telDiam=config.getVal("telDiam")

            self.poke=0
            self.poking=0

            self.control={"poke":0,"close_dm":1,"zero_dm":0,"noiseCovariance":0,"phaseCovariance":0,
                          "takeRef":0,"subtractRef":1}
            self.refCentroids=self.config.getVal("refCentroids",default=None,raiseerror=0)
                                                   # (can be an array or a filename or None.)
            self.takingRef=0
            if type(self.refCentroids)==type(""):
                print("INFORMATION:**tomoRecon**:Using reference centroids "+
                     "from :"+self.refCentroids)
                import util.FITS
                try:
                    self.refCentroids=util.FITS.Read(self.refCentroids)[1]
                except:
##(old)                    print("**************************************************")
                    print("WARNING:**tomoRecon**:**Unable to load reference "+
                        "centroids**: "+self.refCentroids)
##(old)                    print("**************************************************")
                    self.refCentroids=None
            self.saveRefCentroids=self.config.getVal("saveRefCentroids",default=None,raiseerror=0)#either a fits filename or None.
                
            self.phsphs=None
            self.storedPhs=None
            self.montePhaseCovMat=None
            self.subapVariance={}
            self.phssquared=None
            self.lastControlCovariance=0

            self.pokeval=self.config.getVal("pokeval",default=5.)#value to place on actuators when poking.
            self.wfsIDList=[]
            wfsIDList=self.atmosGeom.getWFSOrder(self.idstr[0])#self.config.getVal("wfsIDList",default=
            # self.parent.keys())#the list of parent objects (specifies the order in which the
            # centroids are placed - if that makes a difference?)...
            # print wfsIDList:
            for key in wfsIDList:
                if key not in self.parent.keys():
                    #if get this error, could try specifying a dummy parent that has outputData an array of zeros.
                    raise Exception("tomoRecon: key %s not found in parent, so not using"%str(key))
                else:
                    self.wfsIDList.append(key)
            self.parentDataValid=numpy.zeros((len(self.wfsIDList)),numpy.int32)
            self.multirate=self.config.getVal("multirate",default=0)
            if self.multirate:
                self.dmCommandMulti=numpy.zeros((len(self.wfsIDList),self.nacts),numpy.float32)
            if self.dmModeType=="poke":
                self.nmodes=self.nacts
                self.nLowOrderModalModes=self.config.getVal("nLowOrderModalModes",default=0)#the number
                # of low order modal modes that should be used for each mirror, in addition to the pokes.
                # (this gives a mixed modal/zonal reconstructor).
                if type(self.nLowOrderModalModes)==type(0):
                    self.nLowOrderModalModes=[self.nLowOrderModalModes]*len(self.dmList)
                if len(self.nLowOrderModalModes)!=len(self.dmList):
                    print("ERROR:**tomoRecon**: length of list of number of "+
                        "low order modes for poking doesn't equal number of "+
                        "DMs")
                    self.nLowOrderModalModes=[0]*len(self.dmList)
                n=reduce(lambda x,y:x+y,self.nLowOrderModalModes)
                self.totalLowOrderModalModes=n
                self.nmodes+=n
                self.npokes+=n
                self.modalGain=self.config.getVal("modalGain",default=numpy.ones((n,),numpy.float32))
                self.modeScale=self.config.getVal("modeScale",default=10.)
                self.modalActuatorList=[]
                if n!=0:#modal modes required...
                    modeFile=self.config.getVal("modalFilename",raiseerror=0)
                    writemode="w"
                    #first compute the coords for the actuators within the DMs.
                    for i in range(len(self.dmList)):#for each dm...
                        thisDMnLowOMModes=self.nLowOrderModalModes[i]
                        self.npokesList[i]+=thisDMnLowOMModes # + NAB June/2014
                        dm=self.dmList[i]
                        if dm.zonalDM:
                            dm.computeCoords(2.)
                            coords=dm.coords.copy()
                            coords.shape=(reduce(lambda x,y:x*y,dm.coords.shape[:-1]),dm.coords.shape[-1])
                            coords=numpy.take(coords,numpy.nonzero(self.dmPupList[i].ravel())[0],axis=0) # mod, axis parameter added NAB June/2014
                            nacts=self.nactsList[i]
                            tmp=numpy.zeros((self.nLowOrderModalModes[i],nacts),numpy.float32)
                            self.modalActuatorList.append(tmp)
                            import util.zernikeMod
                            for j in range(0,self.nLowOrderModalModes[i]):#for each zernike...
                                util.zernikeMod.calcZern(j+2,coords,tmp[j])#start at tip/tilt (j+2).
                                tmp[j]*=self.modeScale/numpy.max(numpy.abs(tmp[j]))#so that not too large
                            if modeFile!=None:
                                util.FITS.Write(tmp,modeFile,writeMode=writemode,extraHeader="DMLABEL = %s"%str(dm.label))
                                writemode="a"
                        else:
                            self.modalActuatorList.append(None)
            elif self.dmModeType=="modalPoke":
                self.totalLowOrderModalModes=0
                #poke actuators from a file...
                #This can be used for poking mirror modes or similar...
                self.mirrorModes=self.config.getVal("dmModes")
                if type(self.mirrorModes)==type(""):
                    self.mirrorModes=util.FITS.loadBlockMatrix(self.mirrorModes)
                self.nmodes=self.mirrorModes.shape[0]
                self.npokes=self.nmodes
                #we don't know how the modes are separated out, so just assume an equal number for each dm.
                self.npokesList=[self.npokes//len(self.dmList)]*len(self.dmList)
                self.npokesList[-1]+=self.npokes-sum(self.npokesList)
            self.compressedBits=None#used if the rmx is compressed float format
            self.compressedShape=None
            self.compressedWork=None

            self.outputData=numpy.zeros((self.nacts,),numpy.float64)
            self.reconmxFunction=self.config.getVal("reconmxFunction",default=None,raiseerror=0)#a function that
                               # can be called every iteration to generate a new reconstructor matrix on the fly... 
            self.reconmxFilename=self.config.getVal("reconmxFilename",default=None,raiseerror=0)
            self.gainFactor=self.config.getVal("gainFactor",default=0.5)# used to multiply the estimated phase by.
                                                                        # Can be an array or single value.
            self.decayFactorOpen=self.config.getVal("decayFactorOpen",default=0.)

            self.gains=numpy.zeros((self.nacts,),numpy.float32)
            self.gainList=[]#useful for changing gain of a single dm from the gui
            self.gains[:,]=self.gainFactor
            for i in range(len(self.dmList)):
                dm=self.dmList[i]
                self.gains[self.nactsCumList[i]:self.nactsCumList[i+1]]*=dm.gainAdjustment
                #GUI control only... (makes it easier to change the gain for 1 DM):
                self.gainList.append(self.gains[self.nactsCumList[i]:self.nactsCumList[i+1]])
            # for i in range(len(self.nactsList)):
            #    if not self.closedLoopList[i]:
            #        self.gains[self.nactsCumList[i]:self.nactsCumList[i+1]]=self.gainFactorOpen

            # whether to compute the reconstructor after poking:
            self.computeControl=self.config.getVal("computeControl",default=1)
            self.pmxFilename=self.config.getVal("pmxFilename",raiseerror=0,default=None)
            if self.reconType=="fdpcg":
                self.r0=self.atmosGeom.r0#config.getVal("r0")
                self.l0=self.atmosGeom.l0#config.getVal("l0")
                self.sigList=self.config.getVal("wfsSigList")#list of wfssig for each guide star. todo
                self.monteNoiseCovariance=self.config.getVal("reconMonteCarloNoiseCovariance",default=1)
                self.noiseCov=self.config.getVal("reconNoiseCovariance",default=0.25)
                self.pokemxPCG=None
                self.readNoiseSigma=self.config.getVal("wfs_read_sigma")
                self.readNoiseBg=self.config.getVal("wfs_read_mean")
                self.noiseMatrix=None
                self.phaseCov=None
                print("INFORMATION:**tomoRecon**: creating centroid object "+
                     "for pixel factors")
                import util.centroid#don't know why this is needed - python bug?
                c=None
                if self.monteNoiseCovariance:
                    self.noiseCov=numpy.zeros((self.ncents,),numpy.float32)
                    ns=0
                    for i in range(len(self.ngsList+self.lgsList)):
                        gs=(self.ngsList+self.lgsList)[i]
                        c=util.centroid.centroid(gs.nsubx,self.pupil.fn,readnoise=self.readNoiseSigma,
                                                 addPoisson=1,noiseFloor=self.readNoiseSigma*3+self.readNoiseBg,
                                                 readbg=self.readNoiseBg,sig=self.sigList[i])
                        tmp=c.computeNoiseCovariance(20,convertToRad=1)/c.convFactor**2
                        self.noiseCov[ns:ns+self.ncentList[i]]=numpy.take(tmp.ravel(),
                                                                          self.centIndex[ns:ns+self.ncentList[i]])
                        self.noiseCov[ns+self.ncents/2:ns+self.ncents/2+self.ncentList[i]]=\
                            numpy.take(tmp.ravel(),self.centIndex[ns+self.ncents/2:ns+self.ncents/2+self.ncentList[i]])
                        ns+=self.ncentList[i]
                if c==None:
                    c=util.centroid.centroid(self.ngsList[0].nsubx,self.pupil.fn)
                self.radToPxlFactor=1/c.convFactor**2
                self.phasecov=numpy.zeros((self.nacts,),numpy.float32)
                ns=0
                for dm in self.dmList:
                    if dm.zonalDM==0:
                        raise Exception("TODO - modal DM in tomoRecon")
                    bccb=util.phaseCovariance.bccb(nact=dm.nact,telDiam=self.telDiam,r0=self.r0,l0=self.l0,expand=0)
                    ibccb=util.phaseCovariance.invertBCCB(bccb,dm.nact)
                    self.phasecov[ns:ns+dm.nact**2]=1./ibccb[0]#the diagonal...
                    ns+=dm.nact**2
                self.phasecov*=self.radToPxlFactor#to get into pixels...
                self.pcgOptions={}#self.config.getVal("pcgOptions",default={})
                nlayers=len(self.atmosGeom.layerDict.keys())
                self.convergenceWarning=self.config.getVal("fdpcg_convergenceWarning",default=1)#if 1 will print a warning when convergence not reached.
                self.pcgFakePokeMxVal=self.config.getVal("pcgFakePokeMxVal",default=0.31)#will use an artificial poke matrix if >0.  Probably best to do so.  The value here should be approximately the value that would be given by a real poke matrix, and can have a strong influence on strehl.
                if self.pcgFakePokeMxVal==0.31:
                    print("WARNING:**tomoRecon**: you are using the default "+
                          "value of 0.31 for pcgFakePokeMxVal - this could "+
                          "mean your Strehl is less than it could be.  To "+
                          "get a better value, create a poke matrix and look "+
                          "at the typical values of poked centroids, and use "+
                          "this (and try varying a bit, eg maybe slightly  "+
                          "less...).")
                self.pcgConvergenceValue=self.config.getVal("pcgConvergenceValue",default=0.1)#set this lower (eg 0.001) for slight improvement in reconstruction.  This value is the maximum allowed difference between the pcg solution from the previous iteration, before it is considered to have converged.
                print("INFORMATION:**tomoRecon**: creating pcg object")
                self.pcg=util.tomofdpcg.pcg(self.noiseCov,self.phasecov,None,dmList=self.dmList,ngsList=self.ngsList,lgsList=self.lgsList,minIter=1,maxIter=100,convergenceValue=self.pcgConvergenceValue,convergenceWarning=self.convergenceWarning,fakePokeMxVal=self.pcgFakePokeMxVal,telDiam=self.telDiam)
            elif self.reconType=="spmx":#a broken type...
                self.pokeIgnore=self.config.getVal("pokeIgnore",default=0.01)#ignore pokes below this value when creating a real poke matrix.
                self.createTheoreticalPmx=self.config.getVal("createTheoreticalPmx",default=0)
                self.pTc=None#self.spmx.dot(self.inputData)
                if self.createTheoreticalPmx:
                    import util.createPokeMx#dont know why this is needed - python bug?
                    self.spmx=util.createPokeMx.sparseTomo(ngsList=self.ngsList,lgsList=self.lgsList,dmList=self.dmList,telDiam=self.telDiam).mx
                    self.pTp=util.spmatrix.dotWithSelfTransposed(self.spmx)#this can take 10-20 minutes if large.
                    print("INFORMATION:tomoRecon: doing LU decomposition")
                    self.LUdecomp=scipy.linsolve.splu(self.pTp)#do LU decomposition
                else:
                    self.spmx=dummyClass()#needs a matvec method that returns 0 and takes one input.
                    self.LUdecomp=dummyClass()#needs a solve method that just returns 0 and takes one input.
            elif self.reconType=="spmxSVD":#SVD poke matrix...
                self.pokeIgnore=self.config.getVal("pokeIgnore",default=0.01)
                self.createTheoreticalPmx=self.config.getVal("createTheoreticalPmx",default=0)
                self.svdMin=self.config.getVal("svdMin",default=(0.,0.,0.))#minimum value allowed in u,vt and reconmx before being stripped out.  These factors are only used if svdSparsityFactors>=1.
                if len(self.svdMin)<3:
                    self.svdMin=self.svdMin*3
                self.svdSparsityFactors=self.config.getVal("svdSparsityFactors",default=(1.,1.,1.))#fraction to which the reconmx should be made sparse...  the 3 numbers here are the fraction of svd_u, fraction of svd_vt and fraction of reconmx that should be kept.  These factors are used before svdMin factors if <1.
                if len(self.svdSparsityFactors)<3:
                    self.svdSparsityFactors=self.svdSparsityFactors*3
                self.minEig=self.config.getVal("minEig",default=0.2)
                if self.createTheoreticalPmx:
                    import util.createPokeMx#dont know why this is needed - python bug?
                    self.spmx=util.createPokeMx.sparseTomo(ngsList=self.ngsList,lgsList=self.lgsList,dmList=self.dmList,telDiam=self.telDiam).mx
                    #we now convert to a full matrix, and do SVD decomposition.
                    self.createSVDControl(self.spmx.todense())
                elif self.reconmxFunction!=None:
                    self.reconmx=self.reconmxFunction()
                elif self.reconmxFilename!=None:
                    self.reconmx=self.loadReconmx(self.reconmxFilename)
                else:
                    self.reconmx=0.
                    self.spmx=dummyClass()
            elif self.reconType=="spmxGI":#generalised inverse poke matrix approach.
                self.pokeIgnore=self.config.getVal("pokeIgnore",default=0.01)
                self.createTheoreticalPmx=self.config.getVal("createTheoreticalPmx",default=0)
                self.minEig=self.config.getVal("minEig",default=0.2)
                self.svdSparsityFactors=self.config.getVal("svdSparsityFactors",default=(1.,1.,1.))#fraction to which the reconmx should be made sparse... in this case, first 2 values are ignored.
                if len(self.svdSparsityFactors)<3:
                    self.svdSparsityFactors=self.svdSparsityFactors*3
                if self.createTheoreticalPmx:
                    import util.createPokeMx#dont know why this is needed - python bug?
                    self.spmx=util.createPokeMx.sparseTomo(ngsList=self.ngsList,lgsList=self.lgsList,dmList=self.dmList,telDiam=self.telDiam).mx
                    #we now convert to a full matrix, and do SVD decomposition.
                    self.reconmx=scipy.linalg.pinv(self.spmx.todense(),self.minEig)
                    self.sparsifyMx()
                else:
                    self.spmx=dummyClass()
                    self.reconmx=0.
            elif self.reconType=="svd":#standard dense SVD - for low order systems only...
                self.svdSparsityFactors=None
                self.minEig=self.config.getVal("minEig",default=0.2)
                self.createTheoreticalPmx=self.config.getVal("createTheoreticalPmx",default=0)
                if self.createTheoreticalPmx:
                    import util.createPokeMx#dont know why this is needed - python bug?
                    self.spmx=util.createPokeMx.sparseTomo(ngsList=self.ngsList,lgsList=self.lgsList,dmList=self.dmList,telDiam=self.telDiam).mx
                elif self.reconmxFunction!=None:
                    self.reconmx=self.reconmxFunction()
                else:
                    self.reconmx=self.loadReconmx(self.reconmxFilename)

                #elif self.reconmxFilename!=None and os.path.exists(self.reconmxFilename):
                #    self.reconmx=self.loadReconmx(self.reconmxFilename)
                #else:
                #    self.reconmx=0.
            elif self.reconType in ["pinv","reg","regBig","regSmall","regularised"]:
                self.rcond=self.config.getVal("rcond",default=0.)
                if self.reconmxFunction!=None:
                    self.reconmx=self.reconmxFunction()
                else:
                    self.reconmx=self.loadReconmx(self.reconmxFilename)
            elif self.reconType=="MAP":#MAP reconstruction...
                self.r0=self.atmosGeom.r0#config.getVal("r0")
                self.l0=self.atmosGeom.l0#config.getVal("l0")
                self.nthreads=self.config.getVal("nthreads",default="all")#usually an integer... or "all"
                if self.nthreads=="all":#use all available CPUs...
                    self.nthreads=self.config.getVal("ncpu")#getCpus()
                    print("INFORMATION:tomoRecon: Using {0:d} "+
                        "threads".format(self.nthreads) )
##(old)                    print("INFORMATION:tomoRecon: Using %d threads"%self.nthreads)
                self.rcond=self.config.getVal("rcond",default=0.)
                self.phaseCov=None
                self.monteNoiseCov=None
                self.montePhaseCov=None
                if self.reconmxFunction!=None:
                    self.reconmx=self.reconmxFunction()
                else:
                    self.reconmx=self.loadReconmx(self.reconmxFilename)
                #if self.reconmxFilename!=None and os.path.exists(self.reconmxFilename):
                #    self.reconmx=self.loadReconmx(self.reconmxFilename)
                #else:
                #    self.reconmx=0.
                self.phaseCovFilename=self.config.getVal("phaseCovFilename",raiseerror=0)
                self.noiseCovFilename=self.config.getVal("noiseCovFilename",raiseerror=0)
                self.computePhaseCov=self.config.getVal("computePhaseCov",default=0)
                if self.phaseCovFilename!=None and os.path.exists(self.phaseCovFilename):
                    self.phaseCov=util.FITS.loadBlockMatrix(self.phaseCovFilename)
                    print(("WARNING:tomoRecon: Loading phaseCovariance from "+
                         "{0:s} - check that this is actually what you want..."+
                         "").format(self.phaseCovFilename))
##(old)                         "%s - check that this is actually what you want..."+
##(old)                         "")%self.phaseCovFilename)
                    if self.phaseCov.shape!=(self.nacts,self.nacts):
                        print("WARNING:**tomoRecon**:- tomoRecon - phaseCov shape=%s, not equal to nacts (%d)"%(str(self.phaseCov.shape),self.nacts))
                    h=util.FITS.ReadHeader(self.phaseCovFilename)["parsed"]
                    if h.has_key("r0") and float(h["r0"])!=self.r0:
                        print("INFORMATION:**tomoRecon**:phase cov file is for wrong r0: not using")
                        self.phaseCov=None
                    if h.has_key("l0") and float(h["l0"])!=self.l0:
                        print("INFORMATION:**tomoRecon**:phase cov file is for wrong l0: not using")
                        self.phaseCov=None
                    if h.has_key("r2") and float(h["r2"])!=self.pupil.r2:
                        print("INFORMATION:**tomoRecon**:phase cov file is for wrong r2: not using")
                        self.phaseCov=None
                if self.computePhaseCov and type(self.phaseCov)==type(None):
                    #compute the phase covariance
                    blockList=[]
                    for dm in self.dmList:
                        blockList.append(dm.computePhaseCovariance(self.atmosGeom,self.pupil.r2,self.r0,self.l0,nthreads=self.nthreads,mirrorSurface=dm.getMirrorSurface(),width=-1))
                    self.phaseCov=util.blockMatrix.BlockMatrix()
                    self.phaseCov.assign(blockList)
                    if self.phaseCovFilename!=None:
                        #now save the phase covariance.
                        util.FITS.saveBlockMatrix(self.phaseCov,self.phaseCovFilename,extraHeader=["r0      = %g"%self.r0,"l0      = %g"%self.l0,"r2      = %g"%self.pupil.r2])
                    
                if self.noiseCovFilename!=None:
                    self.noiseCov=util.FITS.loadBlockMatrix(self.noiseCovFilename)
                    if self.noiseCov.shape!=(self.ncents,) and self.noiseCov.shape!=(self.ncents,self.ncents):
                        print("WARNING:**tomoRecon**:- tomoRecon - noiseCov "+
                              "shape=%s, not equal to ncents (%d)"%(
                                 str(self.noiseCov.shape),self.ncents))
                        if len(self.noiseCov.shape)==2:#just take the diagonal-but this could be error for LGS case
                            self.noiseCov=self.noiseCov.diagonal()
                        nsubx=(self.ngsList+self.lgsList)[0].nsubx
                        if self.ncents==self.noiseCov.shape[0]*self.nwfs:
                            print("INFORMATION:**tomoRecon**:"+
                      "assuming identical noise covariance for each wfs")
                            tmp=numpy.zeros((self.ncents),numpy.float32)
                            for i in xrange(self.nwfs):
                                tmp[i*self.noiseCov.shape[0]:(i+1)*self.noiseCov.shape[0]]=self.noiseCov
                            self.noiseCov=tmp
                        elif self.noiseCov.shape[0]==nsubx*nsubx*2:
                            #first need to extract the used values and then use them.
                            print("INFORMATION:**tomoRecon**:"+
                      "extracting noise covariance values from FITS file")
                            xindices=self.centIndex[:self.ncentList[0]]
                            yindices=self.centIndex[self.ncents/2:self.ncents/2+self.ncentList[0]]
                            xnoisecov=numpy.take(self.noiseCov,xindices)
                            ynoisecov=numpy.take(self.noiseCov,yindices)
                            self.noiseCov=numpy.zeros((self.ncents,),numpy.float32)

                            for i in xrange(self.nwfs):
                                self.noiseCov[i*xnoisecov.shape[0]:(i+1)*xnoisecov.shape[0]]=xnoisecov
                                self.noiseCov[self.ncents/2+i*ynoisecov.shape[0]:self.ncents/2+(i+1)*\
                                                  ynoisecov.shape[0]]=ynoisecov
                else:
                    self.noiseCov=0.

                        
                self.influenceScalarProd=None
                self.computeInfluenceScalarProd=self.config.getVal("computeInfluenceScalarProd",default=0)
                #should we compute the DM influence function scalar product, and include
                #in the reconstructor computation?
                if self.computeInfluenceScalarProd:
                    self.computeInfluenceProducts()#this will also compute mirrorScale.
                else:
                    if self.dmModeType=="poke":
                        self.computeMirrorScale()
                    elif self.dmModeType=="modalPoke":
                        self.computeMirrorScale()
                        
                        
                        #for i in xrange(self.nmodes):
                        #    m=self.getMirrorMode(i)
                        #    self.mirrorScale[i]=numpy.sqrt(numpy.sum(m*m))
                self.sumcent=0#used for monte noise covariance.
            elif self.reconType=="pcg":
                self.pcgIters=0#the number of iters required
                self.pcgRegVal=self.config.getVal("pcgRegVal",1e-5)
                self.pcgIsSparse=self.config.getVal("pcgSparse",0)
                self.pcgAfile=self.config.getVal("pcgA",raiseerror=0)
                self.pcgBfile=self.config.getVal("pcgB",raiseerror=0)
                self.pcgX0=None
                self.pcgTol=self.config.getVal("pcgTol",1e-05)
                self.pcgMaxiter=self.config.getVal("pcgMaxiter",30)
                self.pcgPreconditioner=self.config.getVal("pcgPreconditioner",raiseerror=0)
                if type(self.pcgPreconditioner)==type(""):
                    if os.path.exists(self.pcgPreconditioner):
                        self.pcgPreconditioner=util.FITS.loadSparse(self.pcgPreconditioner)
                    else:
                        print("INFORMATION:**tomoRecon**:"+
                          "%s not found"%self.pcgPreconditioner)
                        self.pcgPreconditioner=None
                if type(self.pcgAfile)==type(""):
                    if os.path.exists(self.pcgAfile):
                        self.pcgA=util.FITS.loadSparse(self.pcgAfile)
                    else:
                        print("INFORMATION:**tomoRecon**:"+
                          "%s not found"%self.pcgAfile)
                        self.pcgA=None
                else:
                    self.pcgA=self.pcgAfile
                if type(self.pcgBfile)==type(""):
                    if os.path.exists(self.pcgBfile):
                        self.pcgB=util.FITS.loadSparse(self.pcgBfile)
                    else:
                        print("INFORMATION:**tomoRecon**:"+
                          "%s not found"%self.pcgBfile)
                        self.pcgB=None
                else:
                    self.pcgB=self.pcgBfile
            elif self.reconType=="dicure": # UB, 2012 Aug 3rd
                import util.dicure
                nsubx_tmp = self.config.getVal("wfs_nsubx")
                subapMap = self.pupil.getSubapFlag(nsubx_tmp, self.minarea) # get the subaperture map
                print("INFORMATION:**tomoRecon**:SUBAPMAP:", subapMap)
                self.dicure=util.dicure.DiCuRe( self.dmPupList[0] ) # provide DM actuator map, dmPupList[0].
                                              # If there are more DMs the method will still do something,
                                              # but the result will not make sense, since dicure works only
                                              # for a single DM. (UB, 2012Aug15)
            elif self.reconType=="fewha":
                #austrian wavelet reconstructor.
                import fewhaUser
                self.reconObj=fewhaUser.Fewha(self)


            if self.scaleWhenPoking and self.mirrorScale==None:
                self.computeMirrorScale()
            
            self.abortAfterPoke=self.config.getVal("abortAfterPoke",default=0)
            self.inputData=numpy.zeros(self.ncents,numpy.float32)
            #if self.reconmxFilename==None or type(self.reconmx)==numpy.ndarray:
            #    self.inputData=numpy.zeros(self.ncents,numpy.float32)
            #    #self.inputData.savespace(1)
            #else:
            #    self.inputData=numpy.zeros(self.ncents,numpy.float32)
    # END of __init__

    def finalInitialisation(self):
        print("INFORMATION:**tomoRecon**:"+
              ("decay factors of %s applied where closed loop list %s is 0 - "+
              "is this what you intended?")%(
                  str(self.decayFactorOpen),str(self.closedLoopList)))

    def computeMirrorScale(self):
        self.mirrorScale=numpy.zeros((self.nmodes,),numpy.float32)
        pos=0
        for dm in self.dmList:
            if dm.zonalDM:
                if self.dmModeType=="poke":
                    #make the mirror modes and mirror scale:
                    dm.makeLocalMirrorModes(self.atmosGeom,self.pupil.r2,mirrorSurface=dm.getMirrorSurface())
                    self.mirrorScale[pos:pos+dm.mirrorScale.shape[0]]=dm.mirrorScale
                    pos+=dm.mirrorScale.shape[0]
                elif self.dmModeType=="modalPoke":
                    print "TODO - make mirror modes in tomoRecon"
            else:
                self.mirrorScale[pos:pos+dm.nact]=1
                pos+=dm.nact

    def newParent(self,parent,idstr=None):
        raise Exception("tomoRecon newParent not yet implemented")
    
    def generateNext(self,msg=None):
        """Data coming in from parents can be of 2 types:
         - centroid data (lgs)
         - tilt sensitivity data (ngs)
        Determine which is which from the name of the parent object
        (the dictionary key).  If cent*, is centroid data, if
        tiltsens*, is tilt sensitivity data.
        """
        if self.generate==1:
            if self.newDataWaiting:#this is always 1(!)
                nin=0
                for key in self.parent.keys():
                    if self.parent[key].dataValid==1:
                        nin+=1
##                    else:
##                        print("INFORMATION:**tomoRecon**:tomoRecon: Waiting for data from wfs, but not valid")
                if nin>0:
                    if nin==len(self.parent.keys()):
                        self.dataValid=1
                    elif self.multirate==1:
                        self.dataValid=1#some WFSs are valid
                    else:
                        print("WARNING:**tomoRecon**: got some data but not "+
                           "all, setting dataValid=0")
                        self.dataValid=0
                else:
                    self.dataValid=0
            if self.dataValid:
                self.getInput()
                if self.control["noiseCovariance"]:
                    self.computeMonteNoiseCovariance()
                else:
                    self.lastControlCovariance=0
                if self.parent.has_key("atmos") and self.control["phaseCovariance"]:
                    self.computeMontePhaseCovariance()
                self.calc()# Calculate DM updates
        else:
            self.dataValid=0


    def getInput(self):
        """put all centroids into the input array.
        Note, this is the same order as the ncentList has been defined in...
        """

        cnt=0
        for i in range(len(self.wfsIDList)):
            key=self.wfsIDList[i]
            ns=self.ncentList[i]
            self.parentDataValid[i]=self.parent[key].dataValid
            if self.parent[key].dataValid==1:
                if len(self.parent[key].outputData.shape)==3:
                    self.inputData[cnt:cnt+ns]=numpy.take(self.parent[key].outputData.ravel(),self.centIndex[cnt:cnt+ns])
                    self.inputData[cnt+self.ncents/2:cnt+self.ncents/2+ns]=numpy.take(self.parent[key].outputData.ravel(),self.centIndex[cnt+self.ncents/2:cnt+self.ncents/2+ns])
                else:#wfscent has already cut out the unused subaps... (fullWFSOutput==1 in param file)
                    self.inputData[cnt:cnt+ns]=self.parent[key].outputData[:,0]
                    self.inputData[cnt+self.ncents/2:cnt+self.ncents/2+ns]=self.parent[key].outputData[:,1]
                if type(self.subtractTipTilt)==type({}):
                    stt=self.subtractTipTilt[key]
                else:
                    stt=self.subtractTipTilt
                # this should be used for LGS sensors:
                if stt==-1 or (stt==1 and (self.control["poke"]==0 and self.poking==0 and
                                           self.control["takeRef"]==0 and self.takingRef==0)):
                    #remove the mean slopes...
                    #print("INFORMATION:**tomoRecon**:Removing mean slopes")
                    self.inputData[cnt:cnt+ns]-=self.inputData[cnt:cnt+ns].mean()
                    self.inputData[cnt+self.ncents/2:cnt+self.ncents/2+ns] -= \
                        self.inputData[cnt+self.ncents/2:cnt+self.ncents/2+ns].mean()
            cnt+=ns

    def calc(self):

        if self.takingRef==1:
            self.takingRef=2

        #flatten mirrors then store reference centroids...(with lgs spots, they might not be zero!).
        if self.control["takeRef"]:
            self.control["takeRef"]=0
            self.control["zero_dm"]=1
            self.control["close_dm"]=0
            self.takingRef=1

        if self.control["poke"]:
            self.pokeStartClock=time.clock()
            self.pokeStartTime=time.time()
            self.control["poke"]=0
            self.poking=1
            self.control["close_dm"]=0
            self.pokingDMNo=0#the current DM being poked
            self.pokingActNo=0#the current actuator(s) being poked
            self.pokingDMNoLast=0#the current DM being read by centroids
            self.pokingActNoLast=0#the current actuator(s) being read by cents
            dm=self.dmList[0]
            if dm.pokeSpacing!=None:
                self.pokeActMap=numpy.zeros((dm.nact,dm.nact),numpy.int32)
            if self.reconType=="fdpcg":
                print("INFORMATION:**tomoRecon**:Tomographic FDPCG reconstruction does not do poking")
                self.poking=0
                self.control["close_dm"]=1
            elif self.reconType in ["svd","MAP","pinv","reg","regBig","regSmall","regularised","pcg"]:
                print("INFORMATION:**tomoRecon**:Creating poke matrix of type %s shape %d %d"%(str(self.reconDtype),self.nmodes,self.ncents))
                self.spmx=numpy.zeros((self.nmodes,self.ncents),self.reconDtype)
            elif self.reconType in ["spmx","spmxSVD","spmxGI"]:
                print("INFORMATION:**tomoRecon**:Creating sparse poke matrix")
                self.spmxOld=self.spmx
                if self.sparsePmxType=="lil":
                    self.spmx=scipy.sparse.lil_matrix((self.nmodes,self.ncents),dtype="f")
                elif self.sparsePmxType=="csc":
                    self.spmxData=numpy.zeros((self.ndata,),numpy.float32)
                    self.spmxRowind=numpy.zeros((self.ndata,),numpy.int32)
                    self.spmxIndptr=numpy.zeros((self.ncents+1,),numpy.int32)
                    self.spmx=cmod.svd.sparseMatrixCreate(self.ndata,self.nmodes,self.ncents,
                                                          self.spmxData,self.spmxRowind,self.spmxIndptr)
                else:
                    if self.pmxFilename!=None:
                        self.pmxfname=self.pmxFilename+".mmap"
                    else:
                        self.pmxfname="pmx.mmap"
                    self.spmxValidCnt=0
                    #mmap=cmod.utils.mmapArray(self.pmxfname,(self.nmodes*self.ncents+2880/4,),"f")
                    #TODO: try this sometime. note, it will change the type of mmap (from numpy.ndarray 
                    # to numpy.memmap)
                    mmap=numpy.memmap(self.pmxfname,numpy.float32,"w+",shape=(self.nmodes*self.ncents+2880/4,))
                    hdr=mmap[:2880/4].view("c")
                    self.spmx=mmap[2880/4:]
                    default=1
                    if self.sparsePmxType=="mmap":
                        default=0
                    self.transposeDensePmx=self.config.getVal("transposeDensePmx",default=default)
                    # 1 seems to be fastest for converting to csc, but this may/may not be the case
                    # for when virtual memory (mmap'd) paging is required - ie for large pmx that won't
                    # fit in memory.  Want it set to 0 if intending to use the mmap'd dense matrix.
                    if self.transposeDensePmx==0:
                        self.spmx.shape=(self.nmodes,self.ncents)
                        #self.spmx=cmod.utils.mmapArray(self.pmxfname,(self.nmodes,self.ncents),"f")
                    else:
                        self.spmx.shape=(self.ncents,self.nmodes)
                        #self.spmx=cmod.utils.mmapArray(self.pmxfname,(self.ncents,self.nmodes),"f")
                    hdr[:]=numpy.array(list(util.FITS.MakeHeader(self.spmx.shape,"f",doByteSwap=0)))
                    #self.spmx.savespace(1)
            else:
                raise Exception("Unknown reconType %s for poking"%self.reconType)
            print("INFORMATION:**tomoRecon**:"+
                          "Will be poking for %d integrations"%(self.npokes+1))
        if self.control["zero_dm"]:
            self.control["zero_dm"]=0
            self.outputData[:,]=0.
            if self.reconObj!=None and hasattr(self.reconObj,"resetDM"):
                self.reconObj.resetDM()
        if self.takingRef==2:#the reference centorids should now be ready...
            self.takingRef=0
            self.refCentroids=self.inputData.copy()
            if type(self.saveRefCentroids)==type(""):
                util.FITS.Write(self.refCentroids,self.saveRefCentroids)
                print("INFORMATION:**tomoRecon**:"+
                      "Saving reference centroids to file %s"%(
                           self.saveRefCentroids) )
            else:
                print("INFORMATION:**tomoRecon**:"+
                      "Got reference centroids (not saving)")
        if self.control["subtractRef"] and type(self.refCentroids)!=type(None):
            self.inputData-=self.refCentroids# if this raises an error, it means you specified a reference filename, but it wasn't found.
        if self.poking>0 and self.poking<=self.npokes:
            #set actuator(s) to be set.
            if self.dmModeType=="poke":
                self.outputData[:,]=0.
                if self.poking<=self.npokes-self.totalLowOrderModalModes:
                    # find out which DM we're poking, and whether we're poking individually 
                    # or several actuators at once:
                    dm=self.dmList[self.pokingDMNo]
                    if dm.pokeSpacing!=None:
                        #poking several actuators at once
                        self.pokeActMap[:]=0
                        for i in range(self.pokingActNo/dm.pokeSpacing,dm.nact,dm.pokeSpacing):
                            numpy.put(self.pokeActMap[i],range(self.pokingActNo%dm.pokeSpacing,
                                                               dm.nact,dm.pokeSpacing),1)
                        self.pokeActMapFifo.append(self.pokeActMap*dm.dmflag)
                        self.outputData[self.nactsCumList[self.pokingDMNo]:self.nactsCumList[self.pokingDMNo+1]]=\
                            numpy.take(self.pokeActMap.ravel(),numpy.nonzero(dm.dmflag.ravel()))[0]*\
                            self.pokeval#pokeval added by agb 090130.
                        if self.scaleWhenPoking:#Added 090625 as a test when trying to get MAP to work.
                            self.outputData[self.nactsCumList[self.pokingDMNo]:\
                                                self.nactsCumList[self.pokingDMNo+1]]/=self.mirrorScale[
                                self.nactsCumList[self.pokingDMNo]:self.nactsCumList[self.pokingDMNo+1]]
                    else:#poking individually.
                        self.outputData[self.nactsCumList[self.pokingDMNo]+self.pokingActNo]=self.pokeval
                        # variance of output should be:-
                        #  1/len(outputData)-1/len(outputData)**2.0
                        opvar=self.outputData.var()/(len(self.outputData)**-1.0-len(self.outputData)**-2.0)
                        if abs(opvar-1.0)>1e-8:
                           print("WARNING tomoRecon [{1:3d}:{2:3d}] output data var={0:12.10f}".format(opvar,self.pokingDMNo,self.pokingActNo))
                        #self.outputData[self.poking-1]=self.pokeval
                        if self.reconType=="MAP":#reconstructor assumes modes scaled to orthonormal. but DM may not be
                            self.outputData[self.nactsCumList[self.pokingDMNo]+self.pokingActNo]/=\
                                self.mirrorScale[self.nactsCumList[self.pokingDMNo]+self.pokingActNo]
                    self.pokingActNo+=1
                    if self.pokingActNo==self.npokesList[self.pokingDMNo]:
                        self.pokingDMNo+=1
                        print("INFORMATION:**tomoRecon**:"+
                              "switching to DM no.="+str(self.pokingDMNo))
                        self.pokingActNo=0
                        if self.pokingDMNo<len(self.dmList):
                            dm=self.dmList[self.pokingDMNo]
                            if dm.pokeSpacing!=None:
                                self.pokeActMap=numpy.zeros((dm.nact,dm.nact),numpy.int32)
                        else:
                            self.pokeActMap=None
                else:
                    print("INFORMATION:**tomoRecon**:Now poking mirror modes")
                    mode=self.poking-(self.npokes-self.totalLowOrderModalModes)-1#starts at 0 for the first mode to be poked.
                    #mode=self.totalLowOrderModalModes-1
                    for i in range(len(self.dmList)):#find out which DM we're poking now and get the actuators for it.
                        #print("INFORMATION:**tomoRecon**:Mode: %d, nlomm: %d, dm: %d"%(mode,self.nLowOrderModalModes[i],i))
                        if mode<self.nLowOrderModalModes[i]:
                            #this mode in this DM...
                            #print self.outputData.shape
                            #print self.nactsCumList
                            #print i
                            #print mode
                            if self.modalActuatorList[i]!=None:
                                try:
                                    self.outputData[self.nactsCumList[i]:self.nactsCumList[i+1]]=self.modalActuatorList[i][mode]
                                except:
                                    print("ERROR:**tomoRecon**: assignment "+
                                           "(l815)")
                            break
                        else:
                            #move onto the next DM...
                            mode-=self.nLowOrderModalModes[i]
            elif self.dmModeType=="modalPoke":
                if self.poking<=self.nmodes:
                    self.outputData[:]=self.mirrorModes[self.poking-1]
                    self.pokingActNo+=1
                    if self.pokingActNo==self.npokesList[self.pokingDMNo]:
                        self.pokingDMNo+=1
                        print("INFORMATION:**tomoRecon**: switching to DM no.="+str(self.pokingDMNo))
                        self.pokingActNo=0
                        

        if self.poking>1 and self.poking<=self.npokes+1:
            #then use the centroid values from previous poke to fill the poke matrix.
            if self.poking<=self.npokes-self.totalLowOrderModalModes+1:
                #find out which DM we've just poked, and whether it was poked individually 
                # or several actuators at once:
                dm=self.dmList[self.pokingDMNoLast]
                if dm.pokeSpacing!=None:
                    #poking several at once
                    pokenumber=None
                    #Decide which actuator has produced which centroid values (if any), and then fill the poke matrix.
                    self.fillPokemx(dm,self.pokingDMNoLast)
                else:#poked 1 at once so all centroids belong to this actuator
                    pokenumber=self.nactsCumList[self.pokingDMNoLast]+self.pokingActNoLast#poking-2
                self.pokingActNoLast+=1
                if self.pokingActNoLast==self.npokesList[self.pokingDMNoLast]:
                    self.pokingDMNoLast+=1
                    self.pokingActNoLast=0

            else:#it was a modal poke...
                #pokenumber=self.poking-2+self.nacts-(self.npokes-self.totalLowOrderModalModes)#I think this is right!

                mode=self.poking-2-(self.npokes-self.totalLowOrderModalModes)
                pokenumber=mode+self.nacts
                #print("INFORMATION:**tomoRecon**:poke number %d"%pokenumber)
            if pokenumber!=None:
                if self.reconType in ["spmx","spmxSVD","spmxGI"]:
                    if self.sparsePmxType in ["lil","csc"]:
                        for i,val in enumerate(self.inputData):
                            if numpy.fabs(val)>self.pokeIgnore:
                                if self.sparsePmxType=="lil":
                                    self.spmx[pokenumber,i]=val/self.pokeval
                                elif self.sparsePmxType=="csc":
                                    cmod.svd.sparseMatrixInsert(self.spmx,pokenumber,i,val/self.pokeval)
                    if self.sparsePmxType not in ["lil","csc"]:
                        if self.transposeDensePmx==0:
                            self.spmx[pokenumber]=self.inputData/self.pokeval
                        else:
                            self.spmx[:,pokenumber]=self.inputData/self.pokeval
                        self.spmxValidCnt+=numpy.sum(numpy.fabs(self.inputData)>self.pokeIgnore)
                elif self.reconType in ["svd","MAP","pinv","reg","regularised","regSmall","regBig","pcg"]:
                    self.spmx[pokenumber]=self.inputData/self.pokeval
        if self.poking==self.npokes+1:#finished poking
            self.outputData[:,]=0.
            if self.reconType in ["spmx","spmxSVD","spmxGI"]:
                if self.sparsePmxType=="lil":
                    pass
                elif self.sparsePmxType=="csc":#now convert to a sparse matrix...
                    spmx=scipy.sparse.csc_matrix((numpy.array(self.spmxData,copy=0),
                                                  numpy.array(self.spmxRowind,copy=0),
                                                  numpy.array(self.spmxIndptr,copy=0)),(self.nmodes,self.ncents))
                    cmod.svd.sparseMatrixFree(self.spmx)
                    self.spmx=spmx
                elif self.sparsePmxType=="dense":#convert to sparse...
                    transpose=self.transposeDensePmx
                    data=numpy.zeros((self.spmxValidCnt*2,),numpy.float32)#the *2 is to avoid the sparsity 
                                                          # failing due to undercounting - not sure why this happens - 
                                                          # maybe todo with precision?  Maybe not...
                    rowind=numpy.zeros((self.spmxValidCnt*2,),numpy.int32)
                    indptr=numpy.zeros((self.ncents+1,),numpy.int32)
                    print("INFORMATION:**tomoRecon**:Sparsifying poke matrix")
                    c1=time.clock()
                    t1=time.time()
                    if cmod.svd.sparsifyGenInv(self.spmx,self.pokeIgnore/self.pokeval,transpose,
                                               self.nmodes,data,rowind,indptr)!=0:
                        raise Exception("Error sparsifying genInv (%gs)"%(time.time()-t1))
                    print(("INFORMATION:**tomoRecon**:- time for sparsifying "+
                           "poke matrix: %gs in clocks or %g seconds")%(
                              (time.clock()-c1),time.time()-t1))
                    spmx=scipy.sparse.csc_matrix((numpy.array(data,copy=0),numpy.array(rowind,copy=0),
                                                  numpy.array(indptr,copy=0)),(self.nmodes,self.ncents))
                    self.spmx=spmx
            elif self.reconType in ["svd","MAP","pinv","reg","regularised","regSmall","regBig","pcg"]:
                print("INFORMATION:**tomoRecon**:Finished poking")
                pass#don't need to do owt.

            if self.pmxFilename!=None:
                #save the poke matrix to a file...
                if self.reconType in ["spmx","spmxSVD","spmxGI"]:
                    if self.sparsePmxType in ["lil","csc","dense"]:
                        print("INFORMATION:**tomoRecon**:Saving poke matrix "+
                        "to file %s (shape %s)"%(self.pmxFilename,"hi"))
                        self.spmx=self.spmx.tocsc()
                        self.savecsc(self.spmx,self.pmxFilename)
                    else:
                        import shutil
                        print("INFORMATION:**tomoRecon**:Copying poke matrix "+ 
                              "from %s to %s"%(self.pmxfname,self.pmxFilename))
                        shutil.copyfile(self.pmxfname,self.pmxFilename)
                        
                elif self.reconType in ["pcg","svd","MAP","pinv","reg","regularised","regBig","regSmall"]:
                    print("INFORMATION:**tomoRecon**:Saving poke matrix to "+
                    "file %s (shape %s)"%(
                        self.pmxFilename,str(self.spmx.shape)))
                    util.FITS.Write(self.spmx,self.pmxFilename,extraHeader="NACTLIST= '%s'"%str(self.nactsList))
                #util.FITS.Write(self.spmx.data,self.pmxFilename,extraHeader="SHAPE = %s"%str(self.spmx.shape))
                #util.FITS.Write(self.spmx.rowind,self.pmxFilename,writeMode="a")
                #util.FITS.Write(self.spmx.indptr,self.pmxFilename,writeMode="a")

            if self.computeControl:
                if self.reconType=="spmx":
                    print("INFORMATION:**tomoRecon**:"+
                          "Computing sparse control matrix")
                    self.spmx=self.spmx.tocsc()
                    self.pTp=util.spmatrix.dotWithSelfTransposed(self.spmx)#this can take 10-20 minutes if large.
                    print("INFORMATION:**tomoRecon**:"+
                          "tomoRecon: doing LU decomposition")
                    self.LUdecomp=scipy.linsolve.splu(self.pTp)#do LU decomposition
                elif self.reconType=="spmxSVD":
                    print("INFORMATION:**tomoRecon**:"+
                          "Computing SVD control matrix")
                    try:
                        tmp=self.spmx.todense()
                    except:
                        tmp=self.spmx
                    self.createSVDControl(tmp)
                    if self.reconmxFilename!=None:
                        print("INFORMATION:**tomoRecon**:"+
                          "Writing reconmx to file %s"%self.reconmxFilename)
                        util.FITS.Write(self.reconmx,self.reconmxFilename)

                elif self.reconType=="spmxGI":
                    print("INFORMATION:**tomoRecon**:"+
                          "Computing generalized inverse control matrix")
                    self.reconmx=numpy.array(scipy.linalg.pinv(self.spmx.todense(),self.minEig))
                    self.sparsifyMx()
                elif self.reconType=="svd":
                    print("INFORMATION:**tomoRecon**:"+
                          "Computing svd control matrix")
                    self.createSVDControl(self.spmx)
                    if self.reconmxFilename!=None:
                        print("INFORMATION:**tomoRecon**:"+
                          "Writing reconmx to file %s"%self.reconmxFilename)
                        util.FITS.Write(self.reconmx,self.reconmxFilename)
                elif self.reconType=="pinv":
                    print "Doing pseudo inverse to compute control matrix... (this may take a while for large simulations - consider setting up an alternative approach)"
                    self.reconmx=numpy.linalg.pinv(self.spmx,self.rcond).T.astype(numpy.float32)
                    if self.dmModeType=="modalPoke":
                        #expand the rmx back
                        if self.reconmxFilename!=None:
                            util.FITS.Write(self.reconmx,self.reconmxFilename[:-5]+"modal.fits")
                        self.reconmx=numpy.dot(self.mirrorModes.T,self.reconmx).astype(numpy.float32)
                        print "reconmx shape is now: %s"%str(self.reconmx.shape)
                    if self.reconmxFilename!=None:
                        print("INFORMATION:**tomoRecon**:"+
                          "Writing reconmx to file %s"%self.reconmxFilename)
                        util.FITS.Write(self.reconmx,self.reconmxFilename)
                elif self.reconType=="pcg":
                    self.pcgB=self.spmx
                    self.pcgA=quick.dot(self.spmx,self.spmx.T)
                    self.pcgA.ravel()[::self.pcgA.shape[0]+1]+=self.pcgRegVal
                    if type(self.pcgAfile)==type(""):
                        util.FITS.Write(self.pcgA,self.pcgAfile)
                    if type(self.pcgBfile)==type(""):
                        util.FITS.Write(self.pcgB,self.pcgBfile)
                elif self.reconType=="regBig":
                    self.reconmx=util.regularisation.invert(self.spmx,self.rcond,large=1).T.astype(numpy.float32)
                    if self.reconmxFilename!=None:
                        print("INFORMATION:**tomoRecon**:"+
                              "Writing reconmx to file %s"%self.reconmxFilename)
                        util.FITS.Write(self.reconmx,self.reconmxFilename)
                elif self.reconType in ["reg","regSmall","regularised"]:
                    self.reconmx=util.regularisation.invert(self.spmx,self.rcond).T.astype(numpy.float32)
                    if self.reconmxFilename!=None:
                        print("INFORMATION:**tomoRecon**:"+
                           "Writing reconmx to file %s"%self.reconmxFilename)
                        util.FITS.Write(self.reconmx,self.reconmxFilename)

                elif self.reconType=="MAP":
                    print("INFORMATION:**tomoRecon**:"+
                        "Computing MAP control matrix")
                    self.createMAPControl(self.spmx)
                else:
                    if self.reconObj!=None and hasattr(self.reconObj,"computeControl"):
                        self.reconObj.computeControl(self.spmx)
                    
            print(("INFORMATION:**tomoRecon**:Poking took: %gs in CPU time, "+
                  "or %g seconds")%(
                        (time.clock()-self.pokeStartClock),
                        (time.time()-self.pokeStartTime) ))
            self.poking=0
            if self.abortAfterPoke:
                print("INFORMATION:**tomoRecon**:Finished poking - aborting "+
                     "simulation")
                self.config.abort()
                #if Scientific.MPI.world.size==1:
                #    Scientific.MPI.world.abort()
                #else:
                #    Scientific.MPI.world.abort(0)#this raises an error if python, or aborts correctly if mpipython -
                                                 # however, result is as desired - the simulation finishes.
        if self.poking>0:
            self.poking+=1
        if self.control["close_dm"]:#do the reconstruction.
            self.calc2()
    ## END of calc()

    def calc2(self):
        """Reconstruct using poke matrix or pcg."""
        if type(self.reconmxFunction)!=type(None):
            #call a function which can change the reconstructor on a per-iteration basis:
            self.reconmx=self.reconmxFunction()
        if self.reconType not in ["pcg","dicure","fewha"] and type(self.reconmx)==type(0.):
            #if self.reconmxFilename!=None:
            #    print("INFORMATION:**tomoRecon**:Attempting to load reconmx %s"%self.reconmxFilename)
            self.reconmx=self.loadReconmx(self.reconmxFilename)
        #nsubx=self.wfs_nsubx
        #wfsdata=self.wfsdata
        data=self.inputData#numpy.zeros(wfsdata*2,numpy.Float)
        if self.multirate==0:
            for i in range(len(self.nactsList)):
                dm=self.dmList[i]
                d=self.outputData[self.nactsCumList[i]:self.nactsCumList[i+1]]
                if dm.closedLoop and (dm.polcMatrix is not None):
                    if type(dm.polcMatrix)==type(""):
                        print "Loading polc matrix %s"%dm.polcMatrix
                        dm.polcMatrix=util.FITS.Read(dm.polcMatrix)[1]
                    #With polc, we need to apply (dI +gMP) to the previous actuators, where d is decay factor, M is the rmx and P is the pmx.  
                    #dot the d+gMP matrix with actuators
                    d[:]=numpy.dot(dm.polcMatrix,d)
                elif type(self.decayFactorOpen)==numpy.ndarray and self.decayFactorOpen.size==self.outputData.size:
                    d*=self.decayFactorOpen[self.nactsCumList[i]:self.nactsCumList[i+1]]
                elif dm.decayFactor is not None:
                    d*=dm.decayFactor
                elif not self.closedLoopList[i]:# (open loop actuators...)
                    d*=self.decayFactorOpen
        else:#multirate WFS... update the specific WFSs. (and for now, no polc)
            for i in range(len(self.wfsIDList)):
                if self.parentDataValid[i]:
                    if type(self.decayFactorOpen)==numpy.ndarray and self.decayFactorOpen.size==self.outputData.size:
                        self.dmCommandMulti[i]*=self.decayFactorOpen
                    else:
                        for j in range(len(self.nactsList)):
                            dm=self.dmList[j]
                            if dm.decayFactor is not None:
                                self.dmCommandMulti[i][self.nactsCumList[j]:self.nactsCumList[j+1]]*=dm.decayFactor
                            elif not self.closedLoopList[j]:#open loop
                                self.dmCommandMulti[i][self.nactsCumList[j]:self.nactsCumList[j+1]]*=self.decayFactorOpen

            
        if self.extraActs>0:
            self.outputData[-self.extraActs:]*=self.extraActsDecay
        if self.reconType=="fdpcg":
            self.pcg.solve(data,usePrevious=0)# for some reason, usePrevious=1 will 
                   # cause it to blow up eventually... have no idea why - maybe the pcg
                   # algorithm is not too good at starting points close to the initial.

            self.outputData[:,]+=-self.gainFactor*self.pcg.x
        elif self.reconType=="spmx":#sparse poke matrix reconstruction...
            self.pTc=self.spmx.matvec(numpy.array(data))#dot(data)#do the A^Tb multiplication
            # Solve for x using the LU decomposition:
            self.outputData[:,]+=-self.gainFactor*self.LUdecomp.solve(self.pTc)
        elif self.reconType in ["spmxSVD","spmxGI"]:#SVD reconstruction...
            if self.compressedBits!=None:#reconmx is in compressed float format.
                tmp=-self.doCompressedDot(data)
                if tmp.dtype!=self.outputData.dtype:
                    tmp=tmp.astype(self.outputData.dtype)
            elif type(self.reconmx)==numpy.ndarray:
                if self.multirate==0:
                    if data.shape[0]==self.reconmx.shape[0]:
                        tmp=-quick.dot(data,self.reconmx).astype(self.outputData.dtype)
                    else:
                        #print("INFORMATION:**tomorRecon**:self.reconmx.shape,"+
                        #      "data.shape"+str(self.reconmx.shape)+str(data.shape))
                        tmp=-quick.dot(self.reconmx,data).astype(self.outputData.dtype)
                else:#multi-rate WFSs
                    cnt=0
                    self.outputData[:]=0
                    for i in range(len(self.wfsIDList)):
                        ns=self.ncentList[i]
                        if self.parentDataValid[i]:
                            tmp=quick.dot(data[cnt:cnt+ns],self.reconmx[cnt:cnt+ns])#do the x
                            tmp+=quick.dot(data[cnt+self.ncents/2:cnt+self.ncents/2+ns],self.reconmx[cnt+self.ncents/2:cnt+self.ncents/2+ns])#and add the ys
                            self.dmCommandMulti[i]-=self.gains*tmp
                        cnt+=ns
                        self.outputData+=self.dmCommandMulti[i]
                    tmp=None
                                
            elif (hasattr(scipy.sparse,"csr") and type(self.reconmx)==scipy.sparse.csr.csr_matrix) or \
            (type(self.reconmx)==types.InstanceType or hasattr(self.reconmx,"__module__")) and \
            self.reconmx.__module__ in ["scipy.sparse.sparse"]:
                #print self.reconmx.shape,data.shape
                if self.multirate!=0:
                    raise Exception("multirate not supported for sparse reconstructors at the moment")
                if self.reconmx.shape[0]==data.shape[0]:
                    tmp=-self.reconmx.transpose().dot(data)
                else:
                    tmp=-self.reconmx.dot(data)
            else:
                print("WARNING:**tomoRecon**:tomoRecon: reconmx type not "+
                      "known %s"%str(type(self.reconmx)))
                tmp=None
            if type(tmp)!=type(None):
                try:
                    self.outputData[:,]+=self.gains*tmp[:self.nacts]
                except:
                    print self.outputData.shape,self.gains.shape,tmp[:self.nacts].shape,tmp.shape,self.nacts
                    raise
                #and now apply modal stuff if there was any...
                dmoffset=0
                dmno=0
                for mode in range(self.nmodes-self.nacts):
                    while mode-dmoffset>=self.nLowOrderModalModes[dmno]:
                        #prepare for the next dm...
                        dmoffset+=self.nLowOrderModalModes[dmno]
                        dmno+=1
                    if self.modalActuatorList[dmno]!=None:
                        self.outputData[self.nactsCumList[dmno]:self.nactsCumList[dmno+1]]+=self.\
                            modalGain[mode]*tmp[mode+self.nacts]*self.modalActuatorList[dmno][mode-dmoffset]
        elif self.reconType in ["svd","pinv","reg","regularised","regBig","regSmall"]:
            if self.compressedBits!=None:#reconmx is in compressed float format.
                if self.multirate!=0:
                    raise Exception("multirate wfs not supported here in tomoRecon")
                #tmp=-(self.gains*self.doCompressedDot(data))
                self.outputData-=(self.gains*self.doCompressedDot(data))
            else:
                dorecon=1
                if type(self.reconmx)!=numpy.ndarray or \
                        self.reconmx.shape!=(self.gains.shape[0],data.shape[0]):
                    print("INFORMATION:**tomoRecon**:Reconstructor shape "+
                          "should be "+
                          "(%d,%d)"%(self.gains.shape[0],data.shape[0]))
                    dorecon=0
                    if type(self.reconmx)==numpy.ndarray:
                        print("                         :But current shape is "
                              "%s"%str(self.reconmx.shape))
                        if self.reconmx.shape==(data.shape[0],self.gains.shape[0]):#small...
                            dorecon=1
                            if self.reconmx.size<1e6:#small
                                print("INFORMATION:**tomoRecon**:Transposing "+
                                      "(for small matrices only)")
                                self.reconmx=self.reconmx.T.copy()
                            else:
                                print("INFORMATION:**tomoRecon**:WARNING:"+
                                   "Transposing, but no longer c-contiguous - "+
                                   "performance may be reduced")
                                self.reconmx=self.reconmx.T
                    #tmp=numpy.zeros(self.outputData.shape,self.outputData.dtype)
                if dorecon:
                    if self.multirate==0:
                        #tmp=-(self.gains*quick.dot(self.reconmx,data))#.astype(self.outputData.dtype)
                        self.outputData-=(self.gains*quick.dot(self.reconmx,data))
                    else:
                        cnt=0
                        self.outputData[:]=0
                        for i in range(len(self.wfsIDList)):
                            ns=self.ncentList[i]
                            if self.parentDataValid[i]:
                                tmp=quick.dot(self.reconmx[:,cnt:cnt+ns],data[cnt:cnt+ns],)#do the x
                                tmp+=quick.dot(self.reconmx[:,cnt+self.ncents/2:cnt+self.ncents/2+ns],data[cnt+self.ncents/2:cnt+self.ncents/2+ns])#and add the ys
                                self.dmCommandMulti[i]-=self.gains*tmp
                            cnt+=ns
                            self.outputData+=self.dmCommandMulti[i]
        elif self.reconType=="pcg":
            if self.multirate!=0:
                raise Exception("multirate wfs not supported here in tomoRecon")
            #print data.shape,self.pcgB.shape
            if self.pcgB.shape[1]==data.shape:
                b=quick.dot(self.pcgB,data)
            else:
                b=self.pcgB.T.dot(data)
            #print b.shape,self.pcgA.shape
            self.pcgX0,self.pcgIters=scipy.sparse.linalg.cg(self.pcgA,b,self.pcgX0,
                                                            self.pcgTol,self.pcgMaxiter,M=self.pcgPreconditioner)
            #print self.gains.shape,self.pcgX0.shape,iters,self.pcgX0
            tmp=-self.gains*self.pcgX0
            self.outputData[:]+=tmp
        elif self.reconType=="MAP":
            if self.multirate!=0:
                raise Exception("multirate wfs not supported here in tomoRecon")
            if self.compressedBits!=None:#reconmx is in compressed float format.
                tmp=-(self.gains*self.doCompressedDot(data))
            else:
                tmp=-(self.gains*quick.dot(self.reconmx,data))#.astype(self.outputData.dtype)
            if tmp.dtype!=self.outputData.dtype:
                tmp=tmp.astype(self.outputData.dtype)
            self.outputData[:,]+=tmp
        elif self.reconType=="dicure":
            if self.multirate!=0:
                raise Exception("multirate wfs not supported here in tomoRecon")
            tmp=self.gains*self.dicure.calc( self.inputData )
            self.outputData+=tmp 
        else:#fewha, others.
            if self.multirate!=0:
                raise Exception("multirate wfs not supported here in tomoRecon")
            if self.reconObj!=None and hasattr(self.reconObj,"reconstruct"):
                tmp,partial=self.reconObj.reconstruct(self.inputData)
                #tmp*=self.gains
                if partial:
                    self.outputData+=tmp*self.gains
                else:#the reconstructor does the gain...
                    self.outputData[:]=tmp
    # END of calc2()


    def fillPokemx(self,dm,dmindx):
        """Here, when we've been poking multiple actuators at once, we need to decide which centroids \
           belong to which actuator, and then place them into the poke matrix.
        """
        #First compute the poked actuator coords for this DM.
        if not hasattr(dm,"coords"):
            dm.computeCoords(self.telDiam)
        actuators=self.pokeActMapFifo.pop(0)
        nonzero=numpy.nonzero(actuators.ravel())[0]
        if nonzero.size==0:
            print("INFORMATION:**tomoRecon**:Got no actuators used for this "+
                  "poke")
            return
        pokedCoordsx=numpy.take(dm.coords.ravel(),nonzero*2)#coordinates of poked ones only.
        pokedCoordsy=numpy.take(dm.coords.ravel(),nonzero*2+1)#coordinates of poked ones only.
        cpos=0

        cnt=0
        gsList=self.ngsList+self.lgsList
        for i in range(len(self.wfsIDList)):#for each wfs...
            gs=gsList[i]
            gs.computeCoords(self.telDiam,dm.height)
            key=self.wfsIDList[i]
            ns=self.ncentList[i]
            for subap in self.centIndex[cnt:cnt+ns]/2:#for each subap of this wfs
                x=subap%gs.nsubx
                y=subap/gs.nsubx
                dist=(pokedCoordsx-gs.coords[y,x,0])**2+(pokedCoordsy-gs.coords[y,x,1])**2
                m=numpy.argmin(dist)#this actuator is closest.
                # find the coords for this actuator, and then the number that it corresponds to in the poke matrix.
                indx=nonzero[m]
                pos=dm.dmflag.ravel()[:indx].sum()+self.nactsCumList[dmindx]
                #if pos==196+65 and (cpos in [145,146]):
                #    print (m,indx,pos,x,y,subap,cpos,dist,nonzero)
                #    print self.nactsCumList[dmindx]
                #    print dm.dmflag.shape
                #    print pokedCoordsx
                #    print pokedCoordsy
                #    print gs.coords[y,x]
                #    #util.FITS.Write(dm.dmflag,"dmflag.fits")
                #    if cpos==146:
                #        util.FITS.Write(actuators,"acts.fits")
                # pos is the actuator number in the outputData for this DM
                # This is where the centroid should go in the poke matrix.
                if self.reconType in ["svd","pinv","MAP","reg","regularised","regBig","regSmall","pcg"]:
                    self.spmx[pos,cpos]=self.inputData[cpos]/self.pokeval
                    self.spmx[pos,cpos+self.ncents/2]=self.inputData[cpos+self.ncents/2]/self.pokeval
                elif self.reconType in ["spmx","spmxSVD","spmxGI"]:
                    if self.sparsePmxType=="lil":
                        val=self.inputData[cpos]
                        if numpy.fabs(val)>self.pokeIgnore:
                            self.spmx[pos,cpos]=val/self.pokeval
                        val=self.inputData[cpos+self.ncents/2]
                        if numpy.fabs(val)>self.pokeIgnore:
                            self.spmx[pos,cpos+self.ncents/2]=val/self.pokeval
                    elif self.sparsePmxType=="csc":
                        val=self.inputData[cpos]
                        if numpy.fabs(val)>self.pokeIgnore:
                            cmod.svd.sparseMatrixInsert(self.spmx,pos,cpos,val/self.pokeval)
                        val=self.inputData[cpos+self.ncents/2]
                        if numpy.fabs(val)>self.pokeIgnore:
                            cmod.svd.sparseMatrixInsert(self.spmx,pos,cpos+self.ncents/2,val/self.pokeval)
                    else:
                        val=self.inputData[cpos]
                        val2=self.inputData[cpos+self.ncents/2]
                        if numpy.fabs(val)>self.pokeIgnore:
                            if self.transposeDensePmx==0:
                                self.spmx[pos,cpos]=val/self.pokeval
                            else:
                                self.spmx[cpos,pos]=val/self.pokeval
                            self.spmxValidCnt+=1
                        if numpy.fabs(val2)>self.pokeIgnore:
                            if self.transposeDensePmx==0:
                                self.spmx[pos,cpos+self.ncents/2]=val2/self.pokeval
                            else:
                                self.spmx[cpos+self.ncents/2,pos]=val2/self.pokeval
                            self.spmxValidCnt+=1
                cpos+=1
            cnt+=ns

    def doCompressedDot(self,data,rmx=None,bits=None,shape=None,work=None):
        """Here, rmx is a 1D array in compressed float format, with bits bits per element.
        """
        if rmx==None:
            rmx=self.reconmx
        if bits==None:
            bits=self.compressedBits
        if shape==None:
            shape=self.compressedShape
        if work==None:
            work=self.compressedWork
        expMin=self.compressedExpMin
        expMax=self.compressedExpMax
        if shape[0]==data.shape[0]:
            ncents,nacts=shape
            if work==None:
                mem=getMem("MemFree:")
                nrows=min(mem/4/nacts,ncents)
                print("INFORMATION:**tomoRecon**:Compressed rmx decompressing "+
                     "to %d rowsa at a time"%nrows)
                work=numpy.zeros((nrows,nacts),numpy.float32)
                self.compressedWork=work
            nrows=work.shape[0]
            tmp=numpy.zeros((nacts,),numpy.float32)
            r=work.ravel()
            nsteps=(ncents+nrows-1)/nrows
            for i in xrange(nsteps):
                #uncompress into work
                #if (i+1)*nrows>=ncents:#this is the last one...
                if i==nsteps-1:#this is the last one...
                    end=nrows-(nsteps*nrows-ncents)
                else:
                    end=nrows
                #print("INFORMATION:**tomoRecon**:running %d"%i)
                cmod.utils.uncompressFloatArrayAllThreaded(rmx,r[:end*nacts],bits,expMin,expMax,i*nrows*nacts,8)
                #print("INFORMATION:**tomoRecon**:done")
                tmp+=quick.dot(data[i*nrows:i*nrows+end],work[:end])
        else:
            nacts,ncents=shape
            if work==None:
                mem=getMem("MemFree:")
                nrows=min(mem/4/ncents,nacts)
                print("INFORMATION:**tomoRecon**:Compressed rmx decompressing "+
                     "to %d rowsb at a time"%nrows)
                work=numpy.zeros((nrows,ncents),numpy.float32)
                self.compressedWork=work
            nrows=work.shape[0]
            tmp=numpy.empty((nacts,),numpy.float32)
            r=work.ravel()
            nsteps=(nacts+nrows-1)/nrows
            for i in xrange(nsteps):
                if i==nsteps-1:#last one...
                    end=nrows-(nsteps*nrows-nacts)
                else:
                    end=nrows
                #print("INFORMATION:**tomoRecon**:running %d"%i)
                cmod.utils.uncompressFloatArrayAllThreaded(rmx,r[:end*ncents],bits,expMin,expMax,i*nrows*ncents,8)
                #print("INFORMATION:**tomoRecon**:done")
                tmp[i*nrows:i*nrows+end]=quick.dot(work[:end],data)
        return tmp

    def savecsc(self,csc,filename,hdr=None):
        if type(hdr)==type(""):
            hdr=[hdr]
        elif type(hdr)==type(None):
            hdr=[]
        hdr.append("SHAPE   = %s"%str(csc.shape))
        util.FITS.Write(csc.data[:csc.indptr[-1]],filename,extraHeader=hdr)
        util.FITS.Write(csc.rowind[:csc.indptr[-1]],filename,writeMode="a")
        util.FITS.Write(csc.indptr,filename,writeMode="a")
                
    def loadReconmx(self,reconmxFilename):
        if reconmxFilename==None:
            print("WARNING:**tomoRecon**:Reconmxfilename not specified "+
                  "- using 0.")
            return 0.
        print("INFORMATION:**tomoRecon**:tomoRecon: Loading reconstructor "+
            "(%d,%d) from file: %s"%(self.nmodes,self.ncents,reconmxFilename))
        if os.path.exists(reconmxFilename):
            head=util.FITS.ReadHeader(reconmxFilename)["parsed"]
            if head.has_key("COMPBITS"):
                print("INFORMATION:**tomoRecon**:Reconstructor in compressed "+
                     "FP format")
                #its a compressed floating point format rmx...
                reconmx=util.FITS.Read(reconmxFilename)[1]
                self.compressedBits=int(head["COMPBITS"])
                self.compressedExpMin=int(head.get("EXPMIN",1))
                self.compressedExpMax=int(head.get("EXPMAX",255))
                self.compressedShape=eval(head["SHAPE"])
            else:
                # f=util.FITS.Read(reconmxFilename,savespace=1)
                reconmx=util.FITS.loadSparse(reconmxFilename)
        else:
            print("WARNING:**tomoRecon**: unable to load reconmx "+
                  "%s, using 0. instead"%reconmxFilename)
            reconmx=0.
            #f=[0.]
        #if len(f)==2:
        #    reconmx=f[1].astype("f")
        #elif len(f)>2:
        #    if len(f[1].shape)==1:
        #        reconmx=scipy.sparse.csc_matrix((numpy.array(f[1],copy=0),numpy.array(f[3],copy=0),numpy.array(f[5],copy=0)),eval(f[0]["parsed"]["SHAPE"]))
        #    else:#f[3,5,etc] are probably phase covariance/noise cov etc...
        #        reconmx=f[1]
        return reconmx


    def computeInfluenceProducts(self):
        """Here, compute the scalar products of the influence functions.
        Returns a matrix of size nmodes x nmodes.
        Could also use this to compute mirrorScale?  Since mirrorScale
        would be the diagonal?
        This matrix should be block diagonal.  ie different mirrors
        are orthogonal to each other.  Note, that to do things properly,
        the each mirror mode would be a 1D vector of size sum(dmpup*dmpup)
        summed for all mirrors.  Then each mode would change only the surface
        corresponding to its mirror.
        
        """
        #if self.mirrorModes==None:
        #    self.makeMirrorModes()
        n=self.nmodes#self.mirrorModes.shape[0]
        self.influenceScalarProd=numpy.zeros((n,n),numpy.float32)
        dmposi=0
        actposi=0
        for i in xrange(n):
            dmi=self.dmList[dmposi]
            mmi=dmi.mirrorMode[actposi]
            if dmi.mirrorModeCoords!=None:
                mmci=dmi.mirrorModeCoords[actposi]
            else:
                mmci=None
            actposi+=1
            if actposi==dmi.mirrorModes.shape[0]:
                actposi=0
                dmposi+=1

            for j in xrange(i,n):
                dmj=self.dmList[dmposj]
                mmj=dmj.mirrorMode[actposj]
                if dmj.mirrorModeCoords!=None:
                    mmcj=dmj.mirrorModeCoords[actposj]
                else:
                    mmcj=None
                actposj+=1
                if actposj==dmj.mirrorModes.shape[0]:
                    actposj=0
                    dmposj+=1
                #we're now doing the inner loop, computing the influence scalar product...
                if mmci==None:#mirror mode is whole of mirror
                    if mmcj==None:#mirror mode is whole of mirror
                        self.influenceScalarProd[i,j]=numpy.sum(mmi*mmj)
                    else:#mirror mode is only part of mirror
                        self.influenceScalarProd[i,j]=numpy.sum(mmi[mmcj[1]:mmcj[1]+mmj.shape[0],mmcj[0]:mmcj[0]+mmj.shape[1]]*mmj)
                else:#mirror mode is only part of mirror
                    if mmcj==None:#mirror mode is whole of mirror
                        self.influenceScalarProd[i,j]=numpy.sum(mmi*mmj[mmci[1]:mmci[1]+mmi.shape[0],mmci[0]:mmci[0]+mmi.shape[1]])
                    else:#both are local mirror modes: see where they overlap!
                        ysi=yei=ysj=yej=xsi=xsj=xei=xej=None
                        if mmci[0]<mmcj[0]:
                            if mmci[0]+mmi.shape[1]<=mmcj[0]:
                                pass#no overlap
                            else:#some x overlap
                                xsi=mmcj[0]-mmci[0]#mirrormode i starting point
                                xsj=0
                                if mmci[0]-mmi.shape[1]<mmcj[0]+mmj.shape[1]:
                                    xei=mmi.shape[1]#ending point
                                    xej=mmci[0]+mmi.shape[1]-mmcj[0]
                                else:
                                    xej=mmj.shape[1]
                                    xei=mmcj[0]+mmj.shape[1]-mmci[0]
                        else:
                            if mmcj[0]+mmj.shape[1]<mmci[0]:
                                pass # no overlap
                            else:
                                xsi=0
                                xsj=mmci[0]-mmcj[0]
                                if mmcj[0]-mmj.shape[1]<mmci[0]+mmi.shape[1]:
                                    xej=mmj.shape[1]
                                    xei=mmcj[0]+mmj.shape[1]-mmci[0]
                                else:
                                    xei=mmi.shape[1]
                                    xej=mmci[0]+mmi.shape[1]-mmcj[0]
                        #now do the y coords...
                        if mmci[1]<mmcj[1]:
                            if mmci[1]+mmi.shape[0]<=mmcj[1]:
                                pass#no overlap
                            else:#some x overlap
                                ysi=mmcj[1]-mmci[1]#mirrormode i starting point
                                ysj=0
                                if mmci[1]-mmi.shape[0]<mmcj[1]+mmj.shape[0]:
                                    yei=mmi.shape[0]#ending point
                                    yej=mmci[1]+mmi.shape[0]-mmcj[1]
                                else:
                                    yej=mmj.shape[0]
                                    yei=mmcj[1]+mmj.shape[0]-mmci[1]
                        else:
                            if mmcj[1]+mmj.shape[0]<mmci[1]:
                                pass # no overlap
                            else:
                                ysi=0
                                ysj=mmci[1]-mmcj[1]
                                if mmcj[1]-mmj.shape[0]<mmci[1]+mmi.shape[0]:
                                    yej=mmj.shape[0]
                                    yei=mmcj[1]+mmj.shape[0]-mmci[1]
                                else:
                                    yei=mmi.shape[0]
                                    yej=mmci[1]+mmi.shape[0]-mmcj[1]
                        if None not in [xsi,xsj,xei,xej,ysi,ysj,yei,yej]:
                            self.influenceScalarProd[i,j]=nump.sum(mmi[ysi:yei,xsi:xei]*mmj[ysj:yej,xsj:xej])
                #symmetric, so copy...
                self.influenceScalarProd[j,i]=self.influenceScalarProd[i,j]


        #self.mirrorScale=numpy.sqrt(self.influenceScalarProd.diagonal())
                
    def createMAPControl(self,pokemx,phaseCov=None,noiseCov=None,influenceScalarProd=None):
        """Create a MAP control matrix.  Uses the poke matrix, phase
        covariance and noise covariance.  An optional
        influenceScalarProd matrix can be specified - this is the NxN
        matrix of scalar products of the actuator influence functions
        over the aperture, and is a diagonal for orthonormal functions
        (eg zernikes).

        Calculates:
        rmx=C_f^{-1} C_p A^t(AC_pA^t + C_n)^{-1}
        where C_f is influenceScalarProd, C_p is phase covariance, A is the poke matrix and C_n is noise covariance.
        input pokemx.shape is (nmodes,ncents).
        output reconmx is shape nmodes,ncents.
        """
        pokemxT=pokemx
        pokemx=pokemx.transpose()
        phaseIsBlock=0
        noiseIsBlock=0
        if phaseCov==None:
            phaseCov=self.phaseCov
        if noiseCov==None:
            noiseCov=self.noiseCov
        if influenceScalarProd==None:
            influenceScalarProd=self.influenceScalarProd
        if type(influenceScalarProd)==type(0) and influenceScalarProd==0:
            influenceScalarProd=None
        if isinstance(phaseCov,util.blockMatrix.BlockMatrix):
            phaseIsBlock=1
        if isinstance(noiseCov,util.blockMatrix.BlockMatrix):
            noiseIsBlock=1
        if phaseIsBlock:
            self.apa=quick.dot(util.blockMatrix.dot(pokemx,phaseCov),pokemxT)
        else:
            self.apa=quick.dot(quick.dot(pokemx,phaseCov),pokemxT)
        self.iapan=self.apa.copy()
        if type(noiseCov)==type(0) or type(noiseCov)==type(0.) or type(noiseCov)==type(None):
            if noiseCov==0 or noiseCov==None:
                print("INFORMATION:**tomoRecon**:WARNING - tomoRecon."+
                     "createMAPControl not using noise covariance matrix")
                noiseCov=None
            else:#add const to diag
                self.iapan.ravel()[::self.iapan.shape[0]+1]+=noiseCov
        elif len(noiseCov.shape)==1:#we just have the diagonal...
            self.iapan.ravel()[::self.iapan.shape[0]+1]+=noiseCov
        else:#we have the full matrix (which will be diagonal anyway - except if LGS with TT uncertainty).
            if noiseIsBlock:
                util.blockMatrix.add(self.iapan,noiseCov)
            else:
                self.iapan+=noiseCov
        self.iapan=numpy.linalg.pinv(self.iapan,self.rcond)#shape ncents,ncents
        if phaseIsBlock:
            self.pa=phaseCov.dot(pokemxT)
        else:
            self.pa=quick.dot(phaseCov,pokemxT)#shape nmodes,ncents
        if influenceScalarProd!=None:
            self.pa=quick.dot(influenceScalarProd,self.pa)#shape nmodes,ncents
        self.reconmx=quick.dot(self.pa,self.iapan)#shape nmodes,ncents.

        if self.dmModeType=="poke":
            #need to scale the xinetics_dm output to unity...
            #(for MAP, we are assuming that mirror modes (pokes) are
            #scaled s.t. int(mode**2)==1 when actually xinterp_dm is not
            #scaled like this..).
            #we can do this by altering the reconmx here.
            for i in range(self.nmodes):
                self.reconmx[i]/=self.mirrorScale[i]

        if self.reconmxFilename!=None:
            print("INFORMATION:**tomoRecon**:Writing MAP reconmx to file "+
                  "%s"%self.reconmxFilename)
            util.FITS.Write(self.reconmx,self.reconmxFilename,extraHeader="reconmx")
            if phaseIsBlock:
                for b in phaseCov.blockList:
                    util.FITS.Write(b,self.reconmxFilename,extraHeader="phaseCov",writeMode="a")
            else:
                util.FITS.Write(phaseCov,self.reconmxFilename,writeMode="a",extraHeader="phasecov")
            if noiseIsBlock:
                for b in noiseCov.blockList:
                    util.FITS.Write(b,self.reconmxFilename,extraHeader="noiseCov",writeMode="a")
            elif noiseCov!=None:
                util.FITS.Write(noiseCov,self.reconmxFilename,writeMode="a",extraHeader="noisecov")
            util.FITS.Write(pokemxT,self.reconmxFilename,writeMode="a",extraHeader="pokemx")
            if type(self.mirrorScale)!=type(None):
                util.FITS.Write(self.mirrorScale,self.reconmxFilename,writeMode="a",extraHeader="mirScale")
    def computeMonteNoiseCovariance(self):
        """Note, this isn't a proper way of computing noise
        covariance, as it doesn't take into account the spot
        broadening.  Should actually use something like
        science.centCov"""
        print("DEPRECATED:**tomoRecon**:tomoRecon."+
            "computeMonteNoiseCovariance, use science.centCov")
        if self.sumcent==None:
            print("WARNING:**tomoRecon**: computeMonteNoiseCovariance "+
                  "doesn't actually do what it says, since it "+
                  "doesn't take spot broadening into account - "+
                  "try science.centCov instead")
            self.sumcent=numpy.zeros((self.ncents,),numpy.float64)
            self.sum2cent=numpy.zeros((self.ncents,),numpy.float64)
            self.monteNoiseCov=numpy.zeros((self.ncents,),numpy.float64)
            self.noiseCovCount=0
        else:
            self.noiseCovCount+=1
            self.sumcent+=self.inputData
            self.sum2cent+=self.inputData*self.inputData
        #now compute the variance (diagonal only, since covariances are zero).
        self.monteNoiseCov[:]=self.sum2cent/self.noiseCovCount-self.sumcent*self.sumcent/(
            self.noiseCovCount*self.noiseCovCount)

    def getNoiseCovariance(self,name="centroidNoiseCovariance"):
        """This can be called when a science.centCovariance is also being used and has computed covariances.
        Typically, will be called from the GUI.
        """
        self.monteNoiseCov=self.config.postGet(name)
        
    

    def computeMontePhaseCovariance(self):
        """Not sure how this will work!"""
        print("INFORMATION:**tomoRecon**:TODO computeMontePhaseCovariance()")
        pass
    def makeMonteRecon(self):
        """Make a MAP reconstructor using the monte-carlo generated
        phase covariance and noise covariance."""
        if self.montePhaseCov==None:
            phasecov=self.phaseCov
        else:
            phasecov=self.montePhaseCov
            print("INFORMATION:**tomoRecon**:Using monte phase covariance")
        if self.monteNoiseCov==None:
            noisecov=self.noiseCov
        else:
            noisecov=self.monteNoiseCov
            print("INFORMATION:**tomoRecon**:Using monte noise covariance")
        self.createMAPControl(self.pokemx,phasecov,noisecov)

    def createSVDControl(self,pokemx):
        """create a SVD decomposition of pokemx (dense format).
        """
        print("INFORMATION:**tomoRecon**:tomoRecon.createSVDControl - doing "+
               "svd")
        t1=time.time()
        u,a,vt=scipy.linalg.svd(numpy.array(pokemx))
        self.svd_u=numpy.array(u)
        self.svd_a=numpy.array(a)
        self.svd_vt=numpy.array(vt)
        self.eigenVals=self.svd_a.copy()
        a=self.svd_a
        u=self.svd_u
        vt=self.svd_vt
        n_removed=0
        for i in xrange(a.shape[0]):
            if a[i] < self.minEig:
                a[i]=0.
                n_removed += 1
        t2=time.time()
        print(('INFORMATION:**tomoRecon**:Removed %d modes from control '+
               'matrix (took %g s), ignoring) values below %g')%(
                  n_removed,t2-t1,self.minEig))
        print('INFORMATION:**tomoRecon**:Eigenvalues: '+str(a))
        #now do some sparsifying - for testing purposes... ie remove lowest values...
        if self.svdSparsityFactors!=None:
            if self.svdSparsityFactors[0]>=1.:
                if self.svdMin[0]>0.:
                    u[:,]=numpy.where(numpy.fabs(u)<self.svdMin[0],0.,u).astype(u.typecode())
            else:#keep the fraction given by svdSparsityFactors[0].
                self.sparsifyMx(self.svdSparsityFactors[0],u)
            if self.svdSparsityFactors[1]>=1.:
                if self.svdMin[1]>0.:
                    vt[:,]=numpy.where(numpy.fabs(vt)<self.svdMin[1],0.,vt).astype(vt.typecode())
            else:
                self.sparsifyMx(self.svdSparsityFactors[1],vt)
        ut=numpy.transpose(u)
        v=numpy.transpose(vt)
        #id=numpy.identity(len(a))
        #ai=numpy.multiply(a,id)
        #ai=numpy.where(ai != 0, 1/ai, 0)
        print("INFORMATION:**tomoRecon**:tomoRecon.createSVDControl - doing "+
               "matrix multiplies")
        ai=numpy.where(a!=0,1/a,0)
        print("                         :"+str((ut.shape,ai.shape,v.shape)))
        for i in xrange(min(ut.shape[0],ai.shape[0])):
            ut[i]*=ai[i]
        #the prev loop is the same as numpy.matrixmultiply(ai, ut) if ai is the version created by the next 3 (commented out) lines...
        #ai=numpy.identity(vt.shape[0],"d")
        #ai.flat[0:(vt.shape[0]+1)*a.shape[0]:vt.shape[0]+1]=numpy.where(a!=0,1/a,0)
        #ai=ai[:,:ut.shape[0]]
        #print("INFORMATION:**tomoRecon**:tomoRecon:",v.shape,ai.shape,ut.shape,vt.shape,a.shape,u.shape)
        #print v.shape,ut.shape
        if v.shape[0]>ut.shape[0]:#ncents>nmodes...
            self.reconmx = quick.dot(v[:,:ut.shape[0]], ut)#numpy.matrixmultiply(ai, ut))
        else:
            self.reconmx=quick.dot(v,ut[:vt.shape[1]])
        t3=time.time()
        if self.svdSparsityFactors!=None:
            print("INFORMATION:**tomoRecon**:tomoRecon.createSVDControl - "+
                  "doing possible sparsity checks (matrix multiply took "+
                  "%g s)"%(t3-t2))
            #and do some sparsifying (testing)
            if self.svdSparsityFactors[2]>=1.:
                if self.svdMin[2]>0.:
                    self.reconmx[:,]=numpy.where(numpy.fabs(self.reconmx)<self.svdMin[2],0.,self.reconmx).astype(self.reconmx.typecode())
                else:#make the reconmx a certain faction sparse...
                    self.sparsifyMx()
        s1=reduce(lambda x,y:x*y,self.reconmx.shape)
        s2=s1-numpy.sum(self.reconmx.ravel()==0)
        print(("INFORMATION:**tomoRecon**:tomoRecon.reconmx is %g full "+
              "(%d/%d)")%(float(s2)/s1,s2,s1))

    def sparsifyMx(self,frac=None,mx=None):
        """Sparsity (set to zero) a certain fraction of the reconmx."""
        if frac==None:
            frac=self.svdSparsityFactors[2]
        if frac<0. or frac>=1.:
            return
        if type(mx)==type(None):
            mx=self.reconmx
        s1=reduce(lambda x,y:x*y,mx.shape)
        n=int((1-frac)*s1)#the index which is allowed.
        val=numpy.sort(numpy.fabs(mx.ravel()))[n]
        mx[:,]=numpy.where(numpy.fabs(mx)<val,0.,mx).astype(mx.typecode())
        s2=s1-numpy.sum(mx.ravel()==0)
        print(("INFORMATION:**tomoRecon**:tomoRecon.sparsifyMx: Matrix is %g "+
               "full (%d/%d)")%(float(s2)/s1,s2,s1))

        

    def expandPmx(self):
        """expand self.pokemx into self.pokemxPCG which includes all centroids and actuators.
        I think this is broken - needs updating for tomographic..."""
        actpos=0
        for i in range(self.nact):
            for j in range(self.nact):
                if self.dmflag[i,j]==1:
                    self.pokemxPCGtmp[i*self.nact+j]=self.pokemx[actpos]
                    actpos+=1
        cpos=0
        nsa=self.wfs_nsubx*self.wfs_nsubx
        nc=self.ncents/2
        for i in range(self.wfs_nsubx):
            for j in range(self.wfs_nsubx):
                if self.subflag[i,j]==1:
                    self.pokemxPCG[i*self.wfs_nsubx+j,:,]=self.pokemxPCGtmp[:,cpos]
                    self.pokemxPCG[i*self.wfs_nsubx+j+nsa,:,]=self.pokemxPCGtmp[:,cpos+nc]
                    cpos+=1



    def getLayer(self,n):
        dm=self.dmList[n]
        
        data=self.outputData[self.nactsCumList[n]:self.nactsCumList[n+1]]
        if dm.zonalDM:#reshape the data
            tmp=numpy.zeros((dm.nact,dm.nact),numpy.float32)
            indices=numpy.nonzero(self.dmPupList[n].ravel())[0]
            numpy.put(tmp.ravel(),indices,data)
            data=tmp
        else:
            pass
        return data
            
    def displayGSOverlap(self):
        import util.guideStar
        util.guideStar.displayGSOverlap(gsList=self.ngsList+self.lgsList,layerList=[dm.height for dm in self.dmList],telDiam=self.telDiam,telSec=self.config.getVal("telSec"),fill=True,telCol="red",tells="solid",title=0,outfile=None,sourcedir=None,nx=None,scale=2)
    



    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr==None or self.idstr==[None]:
            id=""
        else:
            id=" (%s)"%self.idstr
        txt="""<plot title="%s Output data%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt"/>\n"""%(self.objID,id,objname)
        txt+="""<plot title="%s Use reference centroids%s" ret="data" when="cmd" texttype="1" wintype="mainwindow">\n<cmd>\ndata=%s.control["subtractRef"]\n%s.control["subtractRef"]=1-data\n</cmd>\nbutton=1-data\n</plot>\n"""%(self.objID,id,objname,objname)
        for i in range(len(self.dmList)):#plots of phase for each layer
            dm=self.dmList[i]
            txt+="""<plot title="%s layer%d%s" cmd="data=-%s.getLayer(%d)" ret="data" type="pylab" when="rpt"/>\n"""%(self.objID,i,id,objname,i)

        if self.reconType=="spmx":
            txt+="""<plot title="%s spmx%s" cmd="data=%s.spmx.todense()" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
            txt+="""<plot title="%s spmxTspmx%s" cmd="data=%s.pTp.todense()" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
        elif self.reconType=="spmxSVD":
            txt+="""<plot title="%s SVD evals%s" cmd="data=%s.eigenVals" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
            txt+="""<plot title="%s reconmx%s" cmd="data=%s.reconmx.todense()" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
        elif self.reconType=="svd":
            txt+="""<plot title="%s SVD evals%s" cmd="data=%s.eigenVals" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
            txt+="""<plot title="%s reconmx%s" cmd="data=%s.reconmx" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
            txt+="""<plot title="%s pmx%s" cmd="data=%s.spmx" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
        elif self.reconType in ["pinv","reg","regularised","regBig","regSmall"]:
            txt+="""<plot title="%s reconmx%s" cmd="data=%s.reconmx" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
            txt+="""<plot title="%s pmx%s" cmd="data=%s.spmx" ret="data" type="pylab" when="cmd"/>\n"""%(self.objID,id,objname)
        txt+="""<plot title="%s inputData%s" cmd="data=%s.inputData" ret="data" type="pylab" when="rpt"/>\n"""%(self.objID,id,objname)
        if self.reconObj!=None and hasattr(self.reconObj,"plottable"):
            txt+=self.reconObj.plottable(objname+".reconObj")
        txt+="""<plot title="Display GS overlap" cmd="%s.displayGSOverlap();data='Ensure that running simulation has X forwarding, and close plot before continuing (it freezes simulation)'" type="text" when="cmd"/>\n"""%objname  
            
        return txt
    

def getMem(memtxt="MemTotal:"):
    lines=open("/proc/meminfo").readlines()
    for line in lines:
        if memtxt in line:
            mem=int(line.split()[1])
            multiplier={"kB":1024,"b":1,"B":1,"mB":1024*1024,"MB":1024**2,"GB":1024**3,"gB":1024**3}
            if line.split()[2] in multiplier.keys():
                mem*=multiplier[line.split()[2]]
            else:
                print("INFORMATION:**tomoRecon**:WARNING - multiplier %s not known for memory")
            print("INFORMATION:**tomoRecon**:Total system memory %d bytes"%mem)
            break
    return mem


if __name__=="__main__":
    print("INFORMATION:**tomoRecon**:Not testing tomoRecon")

class dummyClass:
    def solve(self,val):
        return 0
    def matvec(self,val):
        return 0
